"""Core PCR primer-design logic.

This module contains all of the sequence / primer logic with **no GUI
dependency**, so it can be imported, unit-tested and benchmarked headlessly,
driven from the command line (``primer_cli.py``) or wrapped by the Tkinter GUI
(``enhanced_primer_gui.py``).

Accuracy-relevant design notes
------------------------------
* Tm is reported with the *same* salt / dNTP / DNA concentrations that Primer3
  used to design the primer, so the reported Tm is internally consistent
  (see :class:`ThermoParams`).
* In-silico specificity is a true two-orientation amplicon search: every
  exact binding site of both primers is found on the top strand in both the
  forward and reverse sense, and any pair that can form a product within the
  size window is counted (see :func:`in_silico_pcr`).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
import primer3

logger = logging.getLogger("primer_design")

# Characters that are valid, unambiguous DNA bases.
_NON_DNA_RE = re.compile(r"[^ATGC]")
# As above but keeps length (used for Primer3 templates where coords matter).
_NON_DNA_N_RE = re.compile(r"[^ACGTN]")
_SPLIT_RE = re.compile(r"[,\n\r;\t]+")


# --------------------------------------------------------------------------- #
# Parameter containers
# --------------------------------------------------------------------------- #
@dataclass
class ThermoParams:
    """Reaction conditions shared between design and Tm reporting."""

    mv_conc: float = 50.0      # monovalent cation (mM)
    dv_conc: float = 1.5       # divalent cation (mM)
    dntp_conc: float = 0.6     # dNTP (mM)
    dna_conc: float = 50.0     # primer/oligo (nM)


@dataclass
class PrimerParams:
    """User-tunable primer-design constraints."""

    min_size: int = 18
    opt_size: int = 20
    max_size: int = 25
    min_tm: float = 57.0
    opt_tm: float = 60.0
    max_tm: float = 63.0
    min_gc: float = 40.0
    max_gc: float = 60.0
    product_min: int = 100
    product_max: int = 1000
    gc_clamp: int = 1
    max_poly_x: int = 4
    num_return: int = 1
    thermo: ThermoParams = field(default_factory=ThermoParams)

    def validate(self) -> List[str]:
        """Return a list of human-readable problems (empty == OK)."""
        errs: List[str] = []
        if not (self.min_size <= self.opt_size <= self.max_size):
            errs.append("Primer size must satisfy min <= opt <= max.")
        if not (self.min_tm <= self.opt_tm <= self.max_tm):
            errs.append("Tm must satisfy min <= opt <= max.")
        if not (0 <= self.min_gc <= self.max_gc <= 100):
            errs.append("GC%% must satisfy 0 <= min <= max <= 100.")
        if not (0 < self.product_min <= self.product_max):
            errs.append("Product size must satisfy 0 < min <= max.")
        if self.min_size < 1:
            errs.append("Primer min size must be >= 1.")
        return errs

    def to_primer3_global(self) -> dict:
        """Translate into a Primer3 global-args dictionary."""
        return {
            "PRIMER_OPT_SIZE": self.opt_size,
            "PRIMER_MIN_SIZE": self.min_size,
            "PRIMER_MAX_SIZE": self.max_size,
            "PRIMER_OPT_TM": self.opt_tm,
            "PRIMER_MIN_TM": self.min_tm,
            "PRIMER_MAX_TM": self.max_tm,
            "PRIMER_MIN_GC": self.min_gc,
            "PRIMER_MAX_GC": self.max_gc,
            "PRIMER_NUM_RETURN": max(1, self.num_return),
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.product_min, self.product_max]],
            "PRIMER_GC_CLAMP": self.gc_clamp,
            "PRIMER_MAX_POLY_X": self.max_poly_x,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_MAX_SELF_ANY": 8.0,
            "PRIMER_MAX_SELF_END": 3.0,
            "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
            "PRIMER_PAIR_MAX_COMPL_END": 3.0,
            # Keep reported Tm consistent with design conditions.
            "PRIMER_SALT_MONOVALENT": self.thermo.mv_conc,
            "PRIMER_SALT_DIVALENT": self.thermo.dv_conc,
            "PRIMER_DNTP_CONC": self.thermo.dntp_conc,
            "PRIMER_DNA_CONC": self.thermo.dna_conc,
        }


# --------------------------------------------------------------------------- #
# Sequence helpers
# --------------------------------------------------------------------------- #
def clean_sequence(seq: Optional[str]) -> str:
    """Upper-case and strip everything that is not an unambiguous DNA base."""
    if not seq:
        return ""
    return _NON_DNA_RE.sub("", str(seq).upper())


def calc_gc(seq: Optional[str]) -> float:
    """GC%% of the unambiguous bases of ``seq`` (0.0 for empty/invalid)."""
    clean = clean_sequence(seq)
    if not clean:
        return 0.0
    gc = clean.count("G") + clean.count("C")
    return round(100.0 * gc / len(clean), 2)


def calc_tm(seq: Optional[str], thermo: Optional[ThermoParams] = None) -> Optional[float]:
    """Nearest-neighbour Tm via Primer3, or ``None`` if it cannot be computed."""
    clean = clean_sequence(seq)
    if not clean:
        return None
    thermo = thermo or ThermoParams()
    try:
        return round(
            primer3.calc_tm(
                clean,
                mv_conc=thermo.mv_conc,
                dv_conc=thermo.dv_conc,
                dntp_conc=thermo.dntp_conc,
                dna_conc=thermo.dna_conc,
            ),
            2,
        )
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Tm calculation failed for %r: %s", clean, exc)
        return None


def analyze_oligo(seq: str, thermo: Optional[ThermoParams] = None) -> dict:
    """Return secondary-structure metrics for a single oligo.

    ``hairpin_tm`` and ``homodimer_tm`` are the melting temperatures of the
    most stable self-structures; high values (close to the annealing temp)
    indicate a primer prone to forming structure instead of binding template.
    """
    clean = clean_sequence(seq)
    thermo = thermo or ThermoParams()
    out = {"hairpin_tm": None, "homodimer_tm": None}
    if not clean:
        return out
    kw = dict(
        mv_conc=thermo.mv_conc,
        dv_conc=thermo.dv_conc,
        dntp_conc=thermo.dntp_conc,
        dna_conc=thermo.dna_conc,
    )
    try:
        out["hairpin_tm"] = round(primer3.calc_hairpin(clean, **kw).tm, 2)
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("hairpin calc failed: %s", exc)
    try:
        out["homodimer_tm"] = round(primer3.calc_homodimer(clean, **kw).tm, 2)
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("homodimer calc failed: %s", exc)
    return out


def heterodimer_tm(
    seq_a: str, seq_b: str, thermo: Optional[ThermoParams] = None
) -> Optional[float]:
    """Tm of the most stable duplex between two oligos (primer-dimer risk)."""
    a, b = clean_sequence(seq_a), clean_sequence(seq_b)
    if not a or not b:
        return None
    thermo = thermo or ThermoParams()
    try:
        return round(
            primer3.calc_heterodimer(
                a, b,
                mv_conc=thermo.mv_conc,
                dv_conc=thermo.dv_conc,
                dntp_conc=thermo.dntp_conc,
                dna_conc=thermo.dna_conc,
            ).tm,
            2,
        )
    except Exception as exc:  # pragma: no cover - defensive
        logger.debug("heterodimer calc failed: %s", exc)
        return None


# --------------------------------------------------------------------------- #
# I/O: genome, CDS, GFF3
# --------------------------------------------------------------------------- #
def load_genome(fasta_file: str) -> Dict[str, "SeqIO.SeqRecord"]:
    """Load a (multi-)FASTA genome into a ``{id: SeqRecord}`` dict."""
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    if not genome:
        raise ValueError(f"No sequences found in FASTA file: {fasta_file}")
    return genome


def load_cds_sequences(cds_file: str) -> Dict[str, dict]:
    """Load CDS records keyed by lower-cased gene name / locus tag / id."""
    cds_dict: Dict[str, dict] = {}
    for record in SeqIO.parse(cds_file, "fasta"):
        header = record.description
        gene_name = None
        if "gene=" in header:
            gene_name = header.split("gene=")[1].split()[0].strip("[]")
        elif "locus_tag=" in header:
            gene_name = header.split("locus_tag=")[1].split()[0].strip("[]")
        else:
            gene_name = record.id.split("|")[0] if "|" in record.id else record.id
        if gene_name:
            cds_dict[gene_name.lower()] = {
                "sequence": str(record.seq),
                "id": record.id,
                "description": record.description,
            }
    return cds_dict


def parse_gff3_full(gff_file: str) -> List[dict]:
    """Parse all ``gene`` features from a GFF3 file (0-based start, half-open)."""
    genes: List[dict] = []
    with open(gff_file) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attr = fields[:9]
            locus_tag = "unknown"
            gene_name = None
            for token in attr.split(";"):
                token = token.strip()
                if token.startswith("locus_tag="):
                    locus_tag = token.split("=", 1)[1]
                elif token.startswith("gene="):
                    gene_name = token.split("=", 1)[1]
                elif token.startswith("Name=") and gene_name is None:
                    gene_name = token.split("=", 1)[1]
            try:
                start_i, end_i = int(start) - 1, int(end)
            except ValueError:
                logger.warning("Skipping gene with non-integer coords: %s", line.strip())
                continue
            genes.append(
                {
                    "chrom": chrom,
                    "locus_tag": locus_tag,
                    "gene_name": gene_name,
                    "start": start_i,
                    "end": end_i,
                    "strand": 1 if strand == "+" else -1,
                }
            )
    return genes


def parse_target_names(raw: str) -> List[str]:
    """Parse a free-text gene list, ignoring ``#`` comment lines and blanks.

    This fixes the placeholder bug where example comment text was treated as
    real gene filters.
    """
    if not raw:
        return []
    names: List[str] = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for token in _SPLIT_RE.split(line):
            token = token.strip().lower()
            if token and not token.startswith("#"):
                names.append(token)
    # De-duplicate while preserving order.
    seen = set()
    return [n for n in names if not (n in seen or seen.add(n))]


def filter_genes_by_names(
    all_genes: List[dict], target_names: str
) -> Tuple[List[dict], List[str], List[str]]:
    """Filter genes by name / locus tag (case-insensitive).

    Returns ``(filtered_genes, found_names, not_found_names)``.
    An empty filter returns *all* genes.
    """
    names = parse_target_names(target_names)
    if not names:
        return all_genes, [], []

    wanted = set(names)
    filtered: List[dict] = []
    found: set = set()
    for gene in all_genes:
        gene_name_lower = (gene["gene_name"] or "").lower()
        locus_tag_lower = (gene["locus_tag"] or "").lower()
        if gene_name_lower in wanted or locus_tag_lower in wanted:
            filtered.append(gene)
            found.add(gene_name_lower if gene_name_lower in wanted else locus_tag_lower)

    not_found = [n for n in names if n not in found]
    return filtered, sorted(found), not_found


# --------------------------------------------------------------------------- #
# Sequence extraction
# --------------------------------------------------------------------------- #
def get_gene_sequence(genome: Dict, gene: dict) -> str:
    """Return the (strand-corrected) coding sequence string for ``gene``."""
    rec = genome.get(gene["chrom"])
    if rec is None:
        return ""
    seq = rec.seq[gene["start"]:gene["end"]]
    if gene["strand"] == -1:
        seq = seq.reverse_complement()
    return str(seq)


def extract_sequences_to_fasta(
    genome: Dict, gene_list: List[dict], output_file: str, flank_size: int = 100
) -> None:
    """Write 5' flank, 3' flank and gene region for each gene to a FASTA file."""
    with open(output_file, "w") as out_f:
        for gene in gene_list:
            chrom = gene["chrom"]
            gene_name = gene["gene_name"] or gene["locus_tag"]
            locus = gene["locus_tag"]
            strand = gene["strand"]
            start, end = gene["start"], gene["end"]
            rec = genome.get(chrom)
            if rec is None:
                continue

            five_seq = rec.seq[max(0, start - flank_size):start]
            three_seq = rec.seq[end:min(len(rec.seq), end + flank_size)]
            gene_seq = rec.seq[start:end]
            if strand == -1:
                five_seq, three_seq = three_seq.reverse_complement(), five_seq.reverse_complement()
                gene_seq = gene_seq.reverse_complement()

            def write_record(identifier, seq, region_type):
                out_f.write(f"; ---------- {region_type} ----------\n")
                out_f.write(f">{identifier} | gene={gene_name}\n")
                s = str(seq)
                for i in range(0, len(s), 60):
                    out_f.write(s[i:i + 60] + "\n")
                out_f.write("\n")

            write_record(f"{locus}_5prime_flank", five_seq, "5prime_flank")
            write_record(f"{locus}_3prime_flank", three_seq, "3prime_flank")
            write_record(f"{locus}_gene", gene_seq, "gene")


# --------------------------------------------------------------------------- #
# In-silico specificity
# --------------------------------------------------------------------------- #
def _find_all(sub: str, seq: str) -> List[int]:
    """Return every start index of ``sub`` in ``seq`` (overlapping)."""
    positions: List[int] = []
    if not sub:
        return positions
    start = 0
    while True:
        idx = seq.find(sub, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def prepare_genome(genome: Dict) -> Dict[str, str]:
    """Pre-clean a genome once into ``{chrom: cleaned_upper_string}``.

    Cleaning (a regex pass over each contig) is the dominant per-call cost in
    specificity testing, so do it once and reuse the result across every primer
    pair instead of re-cleaning the whole genome for each one.
    """
    prepared: Dict[str, str] = {}
    for chrom, rec in genome.items():
        prepared[chrom] = rec if isinstance(rec, str) else clean_sequence(str(rec.seq))
    return prepared


def _count_mismatches(a: str, b: str, budget: int) -> int:
    """Hamming distance between two equal-length strings, early-exiting once it
    exceeds ``budget`` (returns ``budget + 1`` in that case)."""
    m = 0
    for x, y in zip(a, b):
        if x != y:
            m += 1
            if m > budget:
                return m
    return m


def _primer_binding_sites(
    top: str, primer: str, seed_len: int, max_mm: int
) -> Tuple[List[Tuple[int, int, int]], List[Tuple[int, int, int]]]:
    """Find every forward- and reverse-acting binding site of one primer.

    Returns ``(fwd_sites, rev_sites)`` where each site is
    ``(left_index, primer_len, mismatches)``.

    A primer only extends if its **3′ end** is well matched, so binding is
    *3′-anchored*: the 3′-most ``seed_len`` bases must match exactly (the seed),
    and the full footprint may differ from the primer by at most ``max_mm``
    bases. Seeding on the 3′ end mirrors the biology (5′ mismatches are
    tolerated, 3′ mismatches abolish extension) and keeps the search fast — the
    exact seed is located with ``str.find`` and only the few candidate
    footprints are scored for mismatches.

    With ``seed_len`` falsy or >= the primer length the whole primer must match
    exactly (``max_mm`` is ignored), reproducing the original exact search.
    """
    L = len(primer)
    rc = str(Seq(primer).reverse_complement())
    if not seed_len or seed_len >= L:
        fwd = [(p, L, 0) for p in _find_all(primer, top)]
        rev = [(p, L, 0) for p in _find_all(rc, top)]
        return fwd, rev

    n = len(top)
    fwd_sites: List[Tuple[int, int, int]] = []
    rev_sites: List[Tuple[int, int, int]] = []

    # Forward-acting: the 3′ seed is the last ``seed_len`` bases of the primer,
    # sitting at the right end of the footprint.
    f_seed = primer[L - seed_len:]
    off = L - seed_len
    pos = top.find(f_seed)
    while pos != -1:
        i = pos - off
        if i >= 0 and i + L <= n:
            mm = _count_mismatches(top[i:i + L], primer, max_mm)
            if mm <= max_mm:
                fwd_sites.append((i, L, mm))
        pos = top.find(f_seed, pos + 1)

    # Reverse-acting: the primer's 3′ end maps to the LEFT end of its
    # reverse-complement on the top strand, so the seed anchors there.
    r_seed = rc[:seed_len]
    pos = top.find(r_seed)
    while pos != -1:
        if pos + L <= n:
            mm = _count_mismatches(top[pos:pos + L], rc, max_mm)
            if mm <= max_mm:
                rev_sites.append((pos, L, mm))
        pos = top.find(r_seed, pos + 1)

    return fwd_sites, rev_sites


def in_silico_pcr(
    forward_primer: str,
    reverse_primer: str,
    genome: Dict,
    min_product: int = 50,
    max_product: int = 5000,
    *,
    seed_len: int = 0,
    max_mismatches: int = 0,
) -> dict:
    """Count predicted amplicons for a primer pair across a loaded genome.

    ``genome`` may be a ``{chrom: SeqRecord}`` dict or a pre-cleaned
    ``{chrom: str}`` dict from :func:`prepare_genome` (preferred for batches —
    it avoids re-cleaning the genome for every primer pair).

    Both primers are searched in *both* senses on every contig:

    * forward-acting site = primer found on the top strand (primes rightward),
    * reverse-acting site = reverse-complement of the primer found on the top
      strand (primes leftward).

    Any forward-acting site upstream of a reverse-acting site that yields a
    product within ``[min_product, max_product]`` is an amplicon. This catches
    off-targets the original one-orientation search missed.

    **Mismatch tolerance (accuracy).** With the default ``seed_len=0`` binding
    is exact. Set ``seed_len`` (≈12 is a good default) and ``max_mismatches`` to
    do a 3′-anchored, mismatch-tolerant search instead: a site counts when its
    3′-most ``seed_len`` bases match exactly and the whole primer differs by at
    most ``max_mismatches`` bases. Real primers prime off-targets that carry
    mismatches in their 5′ half, so exact matching *under*-reports
    non-specificity — this models the biology (3′ end must match to extend)
    while keeping 5′ mismatches tolerated.

    Returns ``{"count": int, "amplicons": [(chrom, start, end, size), ...]}``.
    """
    fwd = clean_sequence(forward_primer)
    rev = clean_sequence(reverse_primer)
    if len(fwd) < 10 or len(rev) < 10:
        return {"count": -1, "amplicons": [], "error": "Primer too short / invalid"}
    if not seed_len:
        max_mismatches = 0

    amplicons: List[Tuple[str, int, int, int]] = []
    for chrom, rec in genome.items():
        top = rec if isinstance(rec, str) else clean_sequence(str(rec.seq))
        if not top:
            continue

        # forward-acting: (5'-end index, primer length, mismatches)
        fwd_sites: List[Tuple[int, int, int]] = []
        # reverse-acting: (footprint left index, primer length, mismatches)
        rev_sites: List[Tuple[int, int, int]] = []
        for p in (fwd, rev):
            f, r = _primer_binding_sites(top, p, seed_len, max_mismatches)
            fwd_sites.extend(f)
            rev_sites.extend(r)

        if not fwd_sites or not rev_sites:
            continue

        rev_sites.sort()
        for f_left, _flen, _fmm in fwd_sites:
            for r_left, r_len, _rmm in rev_sites:
                if r_left < f_left:
                    continue
                product = (r_left + r_len) - f_left
                if product < min_product:
                    continue
                if product > max_product:
                    break  # rev_sites sorted; further ones only larger
                amplicons.append((chrom, f_left, r_left + r_len, product))

    return {"count": len(amplicons), "amplicons": amplicons}


def specificity_label(result: dict) -> str:
    """Human-readable summary of :func:`in_silico_pcr` output."""
    count = result.get("count", -1)
    if count < 0:
        return result.get("error", "Invalid")
    if count == 0:
        return "No amplicons found"
    if count == 1:
        return "Specific (1 amplicon)"
    return f"Non-specific ({count} amplicons)"


# --------------------------------------------------------------------------- #
# Primer design
# --------------------------------------------------------------------------- #
# Named regions of a gene template, ordered 5' -> 3'.
PLACEMENT_REGIONS = ("upstream", "internal", "downstream")


def _region_index(name: str) -> int:
    return PLACEMENT_REGIONS.index(name)


def _p3_template(seq: str) -> str:
    """Upper-case template for Primer3, preserving length (non-ACGT -> N).

    Length must be preserved so that region coordinates stay valid; Primer3
    tolerates ``N`` (governed by ``PRIMER_MAX_NS_ACCEPTED``).
    """
    return _NON_DNA_N_RE.sub("N", str(seq).upper())


def _new_result(seq_id: str, placement: str) -> dict:
    return {
        "gene": seq_id,
        "placement": placement,
        "forward": None,
        "reverse": None,
        "fwd_tm": None,
        "fwd_gc": None,
        "rev_tm": None,
        "rev_gc": None,
        "product_size": None,
        "fwd_hairpin_tm": None,
        "rev_hairpin_tm": None,
        "fwd_homodimer_tm": None,
        "rev_homodimer_tm": None,
        "heterodimer_tm": None,
        "status": "No suitable primers",
    }


def design_primers_for_sequence(
    seq_id: str,
    template: str,
    params: PrimerParams,
    *,
    fwd_region: Optional[Tuple[int, int]] = None,
    rev_region: Optional[Tuple[int, int]] = None,
    product_range: Optional[Tuple[int, int]] = None,
    placement: str = "internal",
) -> dict:
    """Design a single best primer pair for one template.

    Optional ``fwd_region`` / ``rev_region`` are ``(start, length)`` windows
    (in template coordinates) constraining where the left / right primer may be
    placed, via Primer3's ``SEQUENCE_PRIMER_PAIR_OK_REGION_LIST``. This is how
    "forward upstream / reverse downstream" style placements are realised.
    ``product_range`` overrides the parameter product-size window for this call
    (needed when a placement spans flanks larger than the default window).

    Returns a result dict with primer sequences, Tm/GC, product size,
    secondary-structure metrics and a status string. Never raises.
    """
    p3_seq = _p3_template(template)
    clean_len = len(clean_sequence(template))
    result = _new_result(seq_id, placement)

    if clean_len < params.min_size * 2:
        result["status"] = f"Template too short ({clean_len} bp)"
        return result
    eff_min_product = (product_range or (params.product_min, params.product_max))[0]
    if len(p3_seq) < eff_min_product:
        result["status"] = (
            f"Template ({len(p3_seq)} bp) shorter than min product "
            f"({eff_min_product} bp)"
        )
        return result

    seq_args = {"SEQUENCE_ID": seq_id, "SEQUENCE_TEMPLATE": p3_seq}
    if fwd_region is not None or rev_region is not None:
        fl = fwd_region if fwd_region is not None else (-1, -1)
        rl = rev_region if rev_region is not None else (-1, -1)
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = [
            [int(fl[0]), int(fl[1]), int(rl[0]), int(rl[1])]
        ]
    global_args = params.to_primer3_global()
    if product_range is not None:
        global_args["PRIMER_PRODUCT_SIZE_RANGE"] = [[int(product_range[0]), int(product_range[1])]]

    try:
        primers = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as exc:
        result["status"] = f"Primer3 error: {str(exc)[:80]}"
        return result

    fwd = primers.get("PRIMER_LEFT_0_SEQUENCE")
    rev = primers.get("PRIMER_RIGHT_0_SEQUENCE")
    if not fwd or not rev:
        explain = primers.get("PRIMER_PAIR_EXPLAIN", "")
        if explain:
            result["status"] = f"No suitable primers ({explain})"
        return result

    thermo = params.thermo
    fwd_struct = analyze_oligo(fwd, thermo)
    rev_struct = analyze_oligo(rev, thermo)
    result.update(
        {
            "forward": fwd,
            "reverse": rev,
            "fwd_tm": calc_tm(fwd, thermo),
            "fwd_gc": calc_gc(fwd),
            "rev_tm": calc_tm(rev, thermo),
            "rev_gc": calc_gc(rev),
            "product_size": primers.get("PRIMER_PAIR_0_PRODUCT_SIZE"),
            "fwd_hairpin_tm": fwd_struct["hairpin_tm"],
            "rev_hairpin_tm": rev_struct["hairpin_tm"],
            "fwd_homodimer_tm": fwd_struct["homodimer_tm"],
            "rev_homodimer_tm": rev_struct["homodimer_tm"],
            "heterodimer_tm": heterodimer_tm(fwd, rev, thermo),
            "status": "OK",
        }
    )
    return result


def build_gene_template(
    genome: Dict, gene: dict, flank_size: int
) -> Optional[Tuple[str, Dict[str, Tuple[int, int]]]]:
    """Assemble ``5'flank + gene + 3'flank`` (in gene orientation) + region map.

    Returns ``(template, {region_name: (start, length)})`` or ``None`` if the
    contig is missing. Region lengths reflect what was actually available
    (flanks are clipped at contig ends).
    """
    rec = genome.get(gene["chrom"])
    if rec is None:
        return None
    start, end, strand = gene["start"], gene["end"], gene["strand"]
    five = rec.seq[max(0, start - flank_size):start]
    three = rec.seq[end:min(len(rec.seq), end + flank_size)]
    gene_seq = rec.seq[start:end]
    if strand == -1:
        five, three = three.reverse_complement(), five.reverse_complement()
        gene_seq = gene_seq.reverse_complement()
    five, gene_seq, three = str(five), str(gene_seq), str(three)
    template = five + gene_seq + three
    regions = {
        "upstream": (0, len(five)),
        "internal": (len(five), len(gene_seq)),
        "downstream": (len(five) + len(gene_seq), len(three)),
    }
    return template, regions


def placement_combos(mode) -> List[Tuple[str, str]]:
    """Resolve a placement ``mode`` into ``[(fwd_region, rev_region), ...]``.

    ``mode`` may be:

    * ``"internal"``  – both primers inside the gene (default),
    * ``"flanking"``  – forward in the upstream flank, reverse in the downstream
      flank (amplifies the whole gene + flanks, e.g. knockout verification),
    * ``"all"``       – every valid permutation where the forward region is not
      3' of the reverse region (6 combinations),
    * a ``(fwd, rev)`` tuple of region names for a custom placement,
    * a list of such tuples.
    """
    if mode == "internal":
        return [("internal", "internal")]
    if mode == "flanking":
        return [("upstream", "downstream")]
    if mode == "all":
        return [(f, r) for f in PLACEMENT_REGIONS for r in PLACEMENT_REGIONS
                if _region_index(f) <= _region_index(r)]
    if isinstance(mode, tuple):
        return [mode]
    if isinstance(mode, list):
        return list(mode)
    raise ValueError(f"Unknown placement mode: {mode!r}")


def design_for_gene(
    genome: Dict, gene: dict, params: PrimerParams, flank_size: int, mode="internal"
) -> List[dict]:
    """Design primers for one gene under one or more placement modes.

    Returns one result dict per requested placement. For ``mode="internal"``
    with no flanks this is a single internal design (back-compatible).
    """
    name = gene["gene_name"] or gene["locus_tag"]
    combos = placement_combos(mode)
    built = build_gene_template(genome, gene, flank_size)
    if built is None:
        r = _new_result(name, combos[0][0] + "->" + combos[0][1])
        r["status"] = "Contig not found"
        return [r]
    template, regions = built

    results: List[dict] = []
    for fwd_name, rev_name in combos:
        label = f"{fwd_name}->{rev_name}"
        if _region_index(fwd_name) > _region_index(rev_name):
            r = _new_result(name, label)
            r["status"] = "Invalid placement (forward 3' of reverse)"
            results.append(r)
            continue
        freg, rreg = regions[fwd_name], regions[rev_name]
        if freg[1] < params.min_size or rreg[1] < params.min_size:
            r = _new_result(name, label)
            r["status"] = (
                f"Region too small ({fwd_name}={freg[1]}bp, {rev_name}={rreg[1]}bp; "
                f"need >= {params.min_size}). Increase flank size?"
            )
            results.append(r)
            continue
        # Feasible product window for this region pair.
        max_prod = min(len(template), (rreg[0] + rreg[1]) - freg[0])
        product_range = (params.product_min, max_prod)
        results.append(
            design_primers_for_sequence(
                name, template, params,
                fwd_region=freg, rev_region=rreg,
                product_range=product_range, placement=label,
            )
        )
    return results


# Column order shared by CSV writers and the GUI.
RESULT_COLUMNS = [
    "Gene Name", "Placement",
    "Forward Primer", "Fwd Tm", "Fwd GC%",
    "Reverse Primer", "Rev Tm", "Rev GC%",
    "Product Length",
    "Fwd Hairpin Tm", "Rev Hairpin Tm",
    "Fwd SelfDimer Tm", "Rev SelfDimer Tm", "Hetero-dimer Tm",
    "Specificity Check", "Status",
]


def result_to_row(r: dict) -> list:
    """Flatten a design result into a CSV row matching ``RESULT_COLUMNS``."""
    def s(v):
        return "N/A" if v is None else v
    return [
        r["gene"], r.get("placement", "internal"),
        s(r["forward"]), s(r["fwd_tm"]), s(r["fwd_gc"]),
        s(r["reverse"]), s(r["rev_tm"]), s(r["rev_gc"]),
        s(r["product_size"]),
        s(r["fwd_hairpin_tm"]), s(r["rev_hairpin_tm"]),
        s(r["fwd_homodimer_tm"]), s(r["rev_homodimer_tm"]), s(r["heterodimer_tm"]),
        r.get("specificity", "Not tested"),
        r["status"],
    ]


def params_summary(params: PrimerParams, extra: str = "") -> str:
    """One-line parameter banner for the top of a CSV."""
    d = asdict(params)
    base = (
        f"Size {d['min_size']}/{d['opt_size']}/{d['max_size']}, "
        f"Tm {d['min_tm']}/{d['opt_tm']}/{d['max_tm']}, "
        f"GC {d['min_gc']}-{d['max_gc']}%, "
        f"Product {d['product_min']}-{d['product_max']}, "
        f"GC-clamp {d['gc_clamp']}, mv={d['thermo']['mv_conc']}mM"
    )
    return f"Parameters: {base}{(', ' + extra) if extra else ''}"
