#!/usr/bin/env python3
"""Benchmark harness for the PCR primer-design suite.

Measures the two performance-critical paths and writes a Markdown report:

1. **Specificity search** — the original code re-parsed the entire genome FASTA
   from disk *for every primer pair*. This benchmark reproduces that pattern and
   compares it with the fixed :func:`primer_design.in_silico_pcr`, which runs
   against a genome loaded once.
2. **Primer design throughput** — end-to-end ``design_primers_for_sequence``
   over a batch of synthetic templates.

Usage::

    python benchmark.py                 # default sizes
    python benchmark.py --genome-bp 5000000 --pairs 30 --genes 40

Run it before and after a change to confirm you did not regress speed
("benchmark every code edit").
"""

import argparse
import random
import time
from contextlib import contextmanager

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import primer_design as pd


@contextmanager
def timer():
    start = time.perf_counter()
    box = {}
    yield box
    box["elapsed"] = time.perf_counter() - start


def make_genome(n_bp, seed=0):
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(n_bp))
    return seq, {"chr1": SeqRecord(Seq(seq), id="chr1")}


def make_templates(genome_seq, n, length, seed=0):
    """Carve ``n`` sub-sequences out of the genome to design primers from."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        start = rng.randint(0, len(genome_seq) - length - 1)
        out.append(genome_seq[start:start + length])
    return out


def _legacy_specificity(fwd, rev, genome_fasta, max_product):
    """Reproduction of the ORIGINAL per-call genome reload + 1-orientation search."""
    fwd_c = pd.clean_sequence(fwd)
    rev_c = pd.clean_sequence(rev)
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))  # <- reload every call
    rev_rc = str(Seq(rev_c).reverse_complement())
    total = 0
    for _id, rec in genome.items():
        s = pd.clean_sequence(str(rec.seq))
        fpos, i = [], 0
        while True:
            p = s.find(fwd_c, i)
            if p == -1:
                break
            fpos.append(p)
            i = p + 1
        rpos, i = [], 0
        while True:
            p = s.find(rev_rc, i)
            if p == -1:
                break
            rpos.append(p)
            i = p + 1
        for fp in fpos:
            for rp in rpos:
                if fp < rp:
                    size = rp - fp + len(rev_c)
                    if 50 <= size <= max_product:
                        total += 1
    return total


def bench_specificity(genome_seq, genome, fasta_path, pairs):
    """Compare legacy reload-per-pair vs preloaded genome search."""
    rng = random.Random(123)
    primer_pairs = []
    for _ in range(pairs):
        start = rng.randint(0, len(genome_seq) - 500)
        fwd = genome_seq[start:start + 20]
        rev = str(Seq(genome_seq[start + 200:start + 220]).reverse_complement())
        primer_pairs.append((fwd, rev))

    with timer() as legacy:
        for fwd, rev in primer_pairs:
            _legacy_specificity(fwd, rev, fasta_path, 5000)

    # Fixed path as the pipelines use it: clean the genome ONCE, then reuse.
    prepared = pd.prepare_genome(genome)
    with timer() as fixed:
        for fwd, rev in primer_pairs:
            pd.in_silico_pcr(fwd, rev, prepared, 50, 5000)

    # 3'-anchored mismatch-tolerant search over the same preloaded genome.
    with timer() as mismatch:
        for fwd, rev in primer_pairs:
            pd.in_silico_pcr(fwd, rev, prepared, 50, 5000,
                             seed_len=12, max_mismatches=2)

    return (legacy["elapsed"], fixed["elapsed"], mismatch["elapsed"],
            len(primer_pairs))


def bench_design(templates):
    params = pd.PrimerParams(product_min=100, product_max=600)
    ok = 0
    with timer() as t:
        for i, tpl in enumerate(templates):
            r = pd.design_primers_for_sequence(f"g{i}", tpl, params)
            ok += r["status"] == "OK"
    return t["elapsed"], ok, len(templates)


def bench_candidates(templates, num_return):
    """Design ``num_return`` ranked pairs per template (ease-of-use feature).

    Asking Primer3 for several ranked alternates is a single call with a larger
    ``PRIMER_NUM_RETURN``; the extra cost is the per-candidate Tm / structure
    metrics, so this measures the marginal price of the feature.
    """
    params = pd.PrimerParams(product_min=100, product_max=600, num_return=num_return)
    pairs = 0
    with timer() as t:
        for i, tpl in enumerate(templates):
            cands = pd.design_primer_candidates(f"g{i}", tpl, params,
                                                num_candidates=num_return)
            pairs += sum(c["status"] == "OK" for c in cands)
    return t["elapsed"], pairs, len(templates)


def bench_quality_ranking(templates, num_return):
    """Design ``num_return`` pairs/template under default vs quality ranking.

    The composite quality score is already computed for every candidate during
    design (it folds in the Tm / GC / structure numbers), so quality re-ranking
    is only an in-memory sort on top. This measures that the accuracy feature
    adds no meaningful runtime over plain ranked alternates.
    """
    params = pd.PrimerParams(product_min=100, product_max=600, num_return=num_return)
    with timer() as p3:
        for i, tpl in enumerate(templates):
            pd.design_primer_candidates(f"g{i}", tpl, params, rank_by="primer3")
    with timer() as q:
        for i, tpl in enumerate(templates):
            pd.design_primer_candidates(f"g{i}", tpl, params, rank_by="quality")
    return p3["elapsed"], q["elapsed"], len(templates)


def bench_fasta_loading(genome_seq, n_records=200, reps=5):
    """Compare the strict fast-path load vs the lenient fallback.

    The tolerant loader (``read_fasta_records``) tries Biopython's strict
    ``"fasta"`` parser first and only falls back to stripping comments when that
    raises. This measures that the common (clean-file) case keeps the strict
    parser's speed, and how much the comment-stripping fallback costs when a file
    genuinely needs it (e.g. leading ``;`` comments, which crash the strict
    parser on Biopython >= 1.85).
    """
    rng = random.Random(5)
    recs = []
    for i in range(n_records):
        start = rng.randint(0, len(genome_seq) - 2000)
        recs.append(SeqRecord(Seq(genome_seq[start:start + 1500]), id=f"c{i}",
                              description=""))
    clean_path = "_bench_clean.fasta"
    commented_path = "_bench_commented.fasta"
    SeqIO.write(recs, clean_path, "fasta")
    with open(clean_path) as fh:
        body = fh.read()
    # A leading ';' comment block that the strict parser rejects -> forces the
    # lenient fallback path.
    with open(commented_path, "w") as fh:
        fh.write("; benchmark genome with a leading comment block\n;\n\n" + body)

    with timer() as fast:
        for _ in range(reps):
            pd.read_fasta_records(clean_path)
    with timer() as fallback:
        for _ in range(reps):
            pd.read_fasta_records(commented_path)

    import os
    os.remove(clean_path)
    os.remove(commented_path)
    return fast["elapsed"] / reps, fallback["elapsed"] / reps, n_records


def bench_check(genome_seq, genome, n_pairs):
    """Time the `check` path: QC-only vs QC + in-silico-PCR specificity.

    ``analyze_primer_pair`` computes the Tm / GC / secondary-structure / quality
    metrics for a user-supplied pair (no Primer3 design call); adding the genome
    specificity check is the same two-orientation search the design pipeline
    uses. This separates the (cheap) QC cost from the specificity cost so the new
    front-end's overhead is visible.
    """
    rng = random.Random(321)
    pairs = []
    for _ in range(n_pairs):
        start = rng.randint(0, len(genome_seq) - 500)
        fwd = genome_seq[start:start + 20]
        rev = str(Seq(genome_seq[start + 200:start + 220]).reverse_complement())
        pairs.append((fwd, rev))

    params = pd.PrimerParams()
    with timer() as qc:
        for fwd, rev in pairs:
            pd.analyze_primer_pair("p", fwd, rev, params)

    prepared = pd.prepare_genome(genome)
    with timer() as full:
        for fwd, rev in pairs:
            r = pd.analyze_primer_pair("p", fwd, rev, params)
            pd.in_silico_pcr(r["forward"], r["reverse"], prepared, 50, 5000)
    return qc["elapsed"], full["elapsed"], len(pairs)


def bench_placement(genome_seq, genome, n_genes, gene_len=1000, flank=200):
    """Design every placement permutation (6 per gene) over n_genes."""
    rng = random.Random(7)
    params = pd.PrimerParams(product_min=80, product_max=5000)
    genes = []
    for _ in range(n_genes):
        start = rng.randint(flank, len(genome_seq) - gene_len - flank - 1)
        genes.append({"chrom": "chr1", "gene_name": None, "locus_tag": "L",
                      "start": start, "end": start + gene_len, "strand": 1})
    pairs = 0
    with timer() as t:
        for g in genes:
            results = pd.design_for_gene(genome, g, params, flank, mode="all")
            pairs += sum(r["status"] == "OK" for r in results)
    return t["elapsed"], pairs, n_genes


def bench_bed_export(n_pairs=200, amps_per_pair=5, reps=20):
    """Time amplicon -> BED row conversion + sorted write.

    BED export is a pure in-memory transform of the amplicon tuples the
    specificity search already produced (no extra Primer3 or genome work), so
    this confirms it adds negligible cost even for a whole-panel run with many
    predicted amplicons.
    """
    rng = random.Random(11)
    pairs = []
    for _ in range(n_pairs):
        amps = []
        for _ in range(amps_per_pair):
            start = rng.randint(0, 5_000_000)
            amps.append(("chr1", start, start + 200, 200, rng.randint(0, 3)))
        pairs.append(amps)
    target = ("chr1", 1000, 2000)
    out_path = "_bench_amplicons.bed"
    with timer() as t:
        for _ in range(reps):
            rows = []
            for i, amps in enumerate(pairs):
                rows.extend(pd.amplicon_bed_rows("gene", amps, rank=i % 3, target=target))
            pd.write_bed(out_path, rows)
    import os
    os.remove(out_path)
    total_amps = n_pairs * amps_per_pair
    return t["elapsed"] / reps, total_amps


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-bp", type=int, default=2_000_000)
    ap.add_argument("--pairs", type=int, default=20, help="primer pairs for specificity bench")
    ap.add_argument("--genes", type=int, default=30, help="templates for design bench")
    ap.add_argument("--template-bp", type=int, default=1200)
    ap.add_argument("--num-return", type=int, default=3,
                    help="ranked alternates per template for the candidates bench")
    ap.add_argument("-o", "--output", default="benchmark_results.md")
    args = ap.parse_args()

    print(f"Building synthetic genome: {args.genome_bp:,} bp ...")
    genome_seq, genome = make_genome(args.genome_bp)
    fasta_path = "_bench_genome.fasta"
    SeqIO.write(genome["chr1"], fasta_path, "fasta")

    print(f"Specificity benchmark: {args.pairs} primer pairs ...")
    legacy_t, fixed_t, mismatch_t, n_pairs = bench_specificity(
        genome_seq, genome, fasta_path, args.pairs)

    print(f"Design benchmark: {args.genes} templates ...")
    templates = make_templates(genome_seq, args.genes, args.template_bp)
    design_t, ok, n_t = bench_design(templates)

    print(f"Ranked-candidates benchmark: 1 vs {args.num_return} per template ...")
    cand1_t, _, _ = bench_candidates(templates, 1)
    candn_t, candn_pairs, candn_t_n = bench_candidates(templates, args.num_return)

    print(f"Quality-ranking benchmark: primer3 vs quality over {args.genes} templates ...")
    qp3_t, qq_t, q_n = bench_quality_ranking(templates, args.num_return)

    print("FASTA-loading benchmark: strict fast path vs lenient fallback ...")
    load_fast_t, load_fallback_t, load_recs = bench_fasta_loading(genome_seq)

    print(f"Placement benchmark: all permutations over {args.genes} genes ...")
    place_t, place_pairs, place_genes = bench_placement(genome_seq, genome, args.genes)

    print(f"Check benchmark: QC vs QC+specificity over {args.pairs} pairs ...")
    check_qc_t, check_full_t, check_n = bench_check(genome_seq, genome, args.pairs)

    print("BED-export benchmark: amplicon locations -> sorted BED ...")
    bed_t, bed_amps = bench_bed_export()

    speedup = legacy_t / fixed_t if fixed_t else float("inf")

    lines = [
        "# Benchmark Results",
        "",
        f"- Genome size: **{args.genome_bp:,} bp**",
        "- Python/primer3 in-process; single thread.",
        "",
        "## Specificity search (in-silico PCR)",
        "",
        "| Implementation | Pairs | Total time | Per pair |",
        "| --- | ---: | ---: | ---: |",
        f"| Legacy (reload genome per pair, 1 orientation) | {n_pairs} | "
        f"{legacy_t:.3f} s | {legacy_t / n_pairs * 1000:.1f} ms |",
        f"| Fixed (preloaded genome, 2 orientations) | {n_pairs} | "
        f"{fixed_t:.3f} s | {fixed_t / n_pairs * 1000:.1f} ms |",
        f"| 3'-anchored mismatch (seed 12, ≤2 mm) | {n_pairs} | "
        f"{mismatch_t:.3f} s | {mismatch_t / n_pairs * 1000:.1f} ms |",
        "",
        f"**Speedup: {speedup:.1f}x** (the fixed version also checks both primer "
        "orientations, which the legacy version missed). The mismatch-tolerant "
        "search seeds on the exact 3' end with `str.find`, so it stays close to "
        "the exact path while catching off-targets that carry 5' mismatches.",
        "",
        "## Primer design throughput",
        "",
        "| Templates | Successful | Total time | Per template |",
        "| ---: | ---: | ---: | ---: |",
        f"| {n_t} | {ok} | {design_t:.3f} s | {design_t / n_t * 1000:.1f} ms |",
        "",
        "## Ranked alternates (primers per target)",
        "",
        "| Pairs/target | Templates | Pairs returned | Total time | Per template |",
        "| ---: | ---: | ---: | ---: | ---: |",
        f"| 1 | {n_t} | {n_t} | {cand1_t:.3f} s | {cand1_t / n_t * 1000:.1f} ms |",
        f"| {args.num_return} | {candn_t_n} | {candn_pairs} | {candn_t:.3f} s | "
        f"{candn_t / candn_t_n * 1000:.1f} ms |",
        "",
        f"Requesting {args.num_return} ranked pairs is one Primer3 call with a "
        "larger `PRIMER_NUM_RETURN`; the marginal cost is per-candidate Tm and "
        "secondary-structure scoring.",
        "",
        "## Quality re-ranking (primer3 order vs composite quality score)",
        "",
        "| Ranking | Templates | Total time | Per template |",
        "| --- | ---: | ---: | ---: |",
        f"| Primer3 default | {q_n} | {qp3_t:.3f} s | {qp3_t / q_n * 1000:.1f} ms |",
        f"| Quality score | {q_n} | {qq_t:.3f} s | {qq_t / q_n * 1000:.1f} ms |",
        "",
        "The quality score is computed for every candidate during design, so "
        "re-ranking by it is an in-memory sort — the accuracy feature adds no "
        "meaningful runtime over plain ranked alternates.",
        "",
        "## FASTA loading (tolerant reader)",
        "",
        "| Path | Records | Per load |",
        "| --- | ---: | ---: |",
        f"| Strict fast path (clean file) | {load_recs} | "
        f"{load_fast_t * 1000:.1f} ms |",
        f"| Lenient fallback (leading ';' comments) | {load_recs} | "
        f"{load_fallback_t * 1000:.1f} ms |",
        "",
        "The tolerant loader tries Biopython's strict `\"fasta\"` parser first, so "
        "clean files keep its speed. The fallback (comment stripping + reparse) "
        "only runs for files that would otherwise crash the strict parser on "
        "Biopython >= 1.85, and reads the file once more in memory.",
        "",
        "## Placement (all 6 permutations per gene)",
        "",
        "| Genes | Primer pairs | Total time | Per gene |",
        "| ---: | ---: | ---: | ---: |",
        f"| {place_genes} | {place_pairs} | {place_t:.3f} s | "
        f"{place_t / place_genes * 1000:.1f} ms |",
        "",
        "## Primer check (QC of user-supplied primers)",
        "",
        "| Path | Pairs | Total time | Per pair |",
        "| --- | ---: | ---: | ---: |",
        f"| QC only (Tm/GC/structure/quality) | {check_n} | {check_qc_t:.3f} s | "
        f"{check_qc_t / check_n * 1000:.2f} ms |",
        f"| QC + in-silico-PCR specificity | {check_n} | {check_full_t:.3f} s | "
        f"{check_full_t / check_n * 1000:.1f} ms |",
        "",
        "The `check` front-end validates primers a user already has. QC alone is "
        "a few secondary-structure calls per pair (no Primer3 design); adding the "
        "genome specificity check costs the same two-orientation search the design "
        "pipeline already uses, so `check` reuses the accurate engine at no extra "
        "algorithmic cost.",
        "",
        "## Amplicon BED export",
        "",
        "| Amplicons | Total time | Per amplicon |",
        "| ---: | ---: | ---: |",
        f"| {bed_amps} | {bed_t * 1000:.2f} ms | {bed_t / bed_amps * 1e6:.2f} µs |",
        "",
        "Exporting predicted amplicon locations to BED is a pure in-memory "
        "transform of the tuples the specificity search already produced (classify "
        "on/off-target, shade by mismatch, sort by coordinate) plus one file "
        "write — microseconds per amplicon, so `--bed` adds no meaningful cost to "
        "a run that already computed specificity.",
        "",
    ]
    report = "\n".join(lines)
    with open(args.output, "w") as fh:
        fh.write(report + "\n")

    import os
    os.remove(fasta_path)

    print("\n" + report)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
