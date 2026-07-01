"""Correctness tests for :mod:`primer_design`.

Run with ``pytest``. These tests use small synthetic sequences with known
answers so they double as regression guards for the accuracy-critical paths
(in-silico PCR orientation handling, the gene-filter comment bug, Tm/GC).
"""

import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import primer_design as pd  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402


# --------------------------------------------------------------------------- #
# Basic sequence helpers
# --------------------------------------------------------------------------- #
def test_clean_sequence():
    assert pd.clean_sequence("atg-cn N") == "ATGC"
    assert pd.clean_sequence("") == ""
    assert pd.clean_sequence(None) == ""


def test_calc_gc():
    assert pd.calc_gc("GGCC") == 100.0
    assert pd.calc_gc("ATAT") == 0.0
    assert pd.calc_gc("ATGC") == 50.0
    assert pd.calc_gc("") == 0.0          # no divide-by-zero
    assert pd.calc_gc("NNNN") == 0.0      # all ambiguous


def test_calc_tm():
    tm = pd.calc_tm("ATGCATGCATGCATGCATGC")
    assert tm is not None and 40 < tm < 80
    assert pd.calc_tm("") is None


# --------------------------------------------------------------------------- #
# Gene-filter comment bug (the placeholder regression)
# --------------------------------------------------------------------------- #
def test_parse_target_names_ignores_comments():
    # The old default placeholder must NOT be parsed as gene names.
    assert pd.parse_target_names("# Examples:\n# sulA, opgH, galU") == []


def test_parse_target_names_real_input():
    assert pd.parse_target_names("sulA, opgH\nrpoS") == ["sula", "opgh", "rpos"]
    # de-dupes, lower-cases, ignores blanks
    assert pd.parse_target_names("SulA,sula,\n\n") == ["sula"]


def test_filter_genes_by_names():
    genes = [
        {"gene_name": "sulA", "locus_tag": "b0958"},
        {"gene_name": "opgH", "locus_tag": "b1049"},
        {"gene_name": None, "locus_tag": "b9999"},
    ]
    filt, found, not_found = pd.filter_genes_by_names(genes, "sulA, b9999, missing")
    names = {g["locus_tag"] for g in filt}
    assert names == {"b0958", "b9999"}
    assert "missing" in not_found

    # empty filter returns everything
    all_g, _, _ = pd.filter_genes_by_names(genes, "")
    assert len(all_g) == 3


# --------------------------------------------------------------------------- #
# Tolerant FASTA loading (Biopython >= 1.85 strict-parser regression)
# --------------------------------------------------------------------------- #
def test_load_genome_tolerates_leading_comments(tmp_path):
    """A FASTA with blank/';' comment lines before the first record must load.

    Biopython >= 1.85 made the default 'fasta' parser reject such files; the
    loader falls back to a lenient parse so real-world genomes still work.
    """
    fa = tmp_path / "g.fasta"
    fa.write_text("; a comment line\n\n>chr1 desc\nACGTACGTACGTACGT\n")
    genome = pd.load_genome(str(fa))
    assert list(genome) == ["chr1"]
    assert str(genome["chr1"].seq) == "ACGTACGTACGTACGT"


def test_load_genome_plain_fasta_still_works(tmp_path):
    """The strict fast path is still used for clean files (no regression)."""
    fa = tmp_path / "g.fasta"
    fa.write_text(">c1\nACGTACGTACGT\n>c2\nTTTTGGGG\n")
    genome = pd.load_genome(str(fa))
    assert set(genome) == {"c1", "c2"}


def test_load_genome_empty_file_clear_error(tmp_path):
    """An empty/headerless FASTA raises the clear 'No sequences' message, not an
    opaque Biopython traceback."""
    fa = tmp_path / "empty.fasta"
    fa.write_text("\r\n")
    with pytest.raises(ValueError, match="No sequences found"):
        pd.load_genome(str(fa))


def test_load_cds_tolerates_leading_and_inline_comments(tmp_path):
    """CDS loading shares the tolerant reader: leading blanks + ';' comments OK."""
    fa = tmp_path / "cds.fasta"
    fa.write_text("\n\n>g1 [gene=abc]\nATGAAACGT\n; mid comment\n"
                  ">g2 [locus_tag=L2]\nATGCCCGGG\n")
    cds = pd.load_cds_sequences(str(fa))
    assert set(cds) == {"abc", "l2"}
    assert cds["abc"]["sequence"] == "ATGAAACGT"


def test_strip_leading_fasta_comments():
    text = "; hdr comment\n\n>a\nACGT\n; inline\n>b\nTTTT\n"
    cleaned = pd._strip_leading_fasta_comments(text)
    assert cleaned == ">a\nACGT\n>b\nTTTT\n"


# --------------------------------------------------------------------------- #
# GFF3 parsing
# --------------------------------------------------------------------------- #
def test_parse_gff3(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1;gene=abc;locus_tag=L1\n"
        "chr1\tsrc\tgene\t300\t400\t.\t-\t.\tID=g2;Name=def;locus_tag=L2\n"
        "chr1\tsrc\tCDS\t100\t200\t.\t+\t.\tID=c1\n"  # ignored
    )
    genes = pd.parse_gff3_full(str(gff))
    assert len(genes) == 2
    assert genes[0]["start"] == 99 and genes[0]["end"] == 200  # 0-based, half-open
    assert genes[0]["strand"] == 1
    assert genes[1]["gene_name"] == "def"   # falls back to Name=
    assert genes[1]["strand"] == -1


# --------------------------------------------------------------------------- #
# In-silico PCR: orientation correctness (the accuracy fix)
# --------------------------------------------------------------------------- #
def _genome_from_str(seq):
    return {"chr1": SeqRecord(Seq(seq), id="chr1")}


def test_in_silico_pcr_single_amplicon():
    fwd = "AAAGGGCTACGTTGCAATCG"          # 20 bp, non-palindromic
    rev = "CCCTTTGCAGTACTTAGGAC"          # 20 bp, non-palindromic
    rev_rc = str(Seq(rev).reverse_complement())
    spacer = "A" * 100
    template = "GG" + fwd + spacer + rev_rc + "GG"
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=50, max_product=500)
    assert res["count"] == 1
    chrom, start, end, size, mm = res["amplicons"][0]
    assert size == len(fwd) + len(spacer) + len(rev)
    assert mm == 0  # exact match -> no mismatches


def test_in_silico_pcr_detects_second_site():
    """Two binding sites for the pair must be reported as non-specific."""
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    amp = fwd + ("A" * 100) + rev_rc
    template = "GG" + amp + ("T" * 300) + amp + "GG"
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=50, max_product=500)
    assert res["count"] == 2


def test_in_silico_pcr_reverse_orientation():
    """An amplicon where the primers' roles are swapped must still be found.

    The original code only checked fwd-on-top + rev-revcomp-on-top and would
    miss this; the fixed two-orientation search must catch it.
    """
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    fwd_rc = str(Seq(fwd).reverse_complement())
    # Layout: reverse primer acts as the left/forward primer here.
    template = "GG" + rev + ("A" * 100) + fwd_rc + "GG"
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=50, max_product=500)
    assert res["count"] >= 1


def test_in_silico_pcr_product_window():
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    template = fwd + ("A" * 100) + rev_rc
    # Window excludes the real product -> no amplicon.
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=500, max_product=1000)
    assert res["count"] == 0


def test_in_silico_pcr_mismatch_offtarget_detected():
    """A 5'-mismatched off-target (3' seed intact) is invisible to exact search
    but must be caught by the 3'-anchored mismatch-tolerant search."""
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    spacer = "A" * 100
    # Off-target forward site: two 5' mismatches, identical 3' seed.
    fwd_5mm = "CC" + fwd[2:]
    primary = fwd + spacer + rev_rc
    offtarget = fwd_5mm + spacer + rev_rc
    template = "GG" + primary + ("T" * 300) + offtarget + "GG"
    genome = _genome_from_str(template)

    exact = pd.in_silico_pcr(fwd, rev, genome, min_product=50, max_product=200)
    assert exact["count"] == 1  # only the perfect site

    mm = pd.in_silico_pcr(fwd, rev, genome, min_product=50, max_product=200,
                          seed_len=12, max_mismatches=2)
    assert mm["count"] == 2  # perfect site + 5'-mismatched off-target


def test_in_silico_pcr_3prime_mismatch_not_extended():
    """A mismatch inside the 3' seed abolishes priming, so the off-target must
    NOT be counted even though it is only one base off overall."""
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    spacer = "A" * 100
    # Off-target with a single mismatch at the 3'-most base (breaks the seed).
    fwd_3mm = fwd[:-1] + ("A" if fwd[-1] != "A" else "C")
    template = ("GG" + fwd + spacer + rev_rc
                + ("T" * 300) + fwd_3mm + spacer + rev_rc + "GG")
    mm = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                          min_product=50, max_product=200,
                          seed_len=12, max_mismatches=3)
    assert mm["count"] == 1  # 3'-broken site rejected despite mismatch budget


def test_in_silico_pcr_seed_len_zero_is_exact():
    """seed_len=0 must reproduce exact-match behaviour bit-for-bit."""
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    template = "GG" + fwd + ("A" * 100) + rev_rc + "GG"
    genome = _genome_from_str(template)
    a = pd.in_silico_pcr(fwd, rev, genome, 50, 500)
    b = pd.in_silico_pcr(fwd, rev, genome, 50, 500, seed_len=0, max_mismatches=5)
    assert a == b


def test_specificity_label():
    assert pd.specificity_label({"count": 0}) == "No amplicons found"
    assert pd.specificity_label({"count": 1}).startswith("Specific")
    assert pd.specificity_label({"count": 3}).startswith("Non-specific")
    assert pd.specificity_label({"count": -1, "error": "bad"}) == "bad"


def test_specificity_label_reports_offtarget_mismatches():
    """A non-specific result surfaces the nearest off-target's mismatch count."""
    res = {"count": 2, "amplicons": [
        ("chr1", 0, 120, 120, 0),   # intended (perfect)
        ("chr1", 400, 520, 120, 2),  # off-target with 2 mismatches
    ]}
    label = pd.specificity_label(res)
    assert "Non-specific (2 amplicons" in label
    assert "2 mismatches" in label
    # Singular grammar for a one-mismatch off-target.
    res["amplicons"][1] = ("chr1", 400, 520, 120, 1)
    assert "1 mismatch)" in pd.specificity_label(res)


def test_specificity_label_singular_amplicon_when_on_target_only():
    """One perfect amplicon is 'Specific', not a mismatch-annotated non-specific."""
    res = {"count": 1, "amplicons": [("chr1", 0, 120, 120, 0)]}
    assert pd.specificity_label(res) == "Specific (1 amplicon)"


# --------------------------------------------------------------------------- #
# BED export of predicted amplicon locations
# --------------------------------------------------------------------------- #
def test_amplicon_bed_rows_classifies_on_and_off_target():
    """The intended locus is labelled ontarget; a distal amplicon offtarget."""
    amplicons = [
        ("chr1", 100, 220, 120, 0),   # overlaps the target -> ontarget
        ("chr1", 5000, 5120, 120, 2),  # elsewhere -> offtarget, 2 mismatches
        ("chr2", 100, 220, 120, 0),   # other contig -> offtarget
    ]
    rows = pd.amplicon_bed_rows("geneA", amplicons, rank=0,
                                target=("chr1", 150, 180))
    assert len(rows) == 3
    # BED coordinates pass straight through (already 0-based, half-open).
    assert rows[0][:3] == ["chr1", 100, 220]
    names = [r[3] for r in rows]
    assert names[0] == "geneA_rank1_ontarget_mm0"
    assert names[1] == "geneA_rank1_offtarget_mm2"
    assert names[2] == "geneA_rank1_offtarget_mm0"  # right contig required
    # Score shades by mismatch: perfect = 1000, each mismatch -250.
    assert rows[0][4] == 1000
    assert rows[1][4] == 500
    assert all(r[5] == "+" for r in rows)


def test_amplicon_bed_rows_no_target_is_all_offtarget():
    rows = pd.amplicon_bed_rows("g", [("c", 0, 100, 100, 0)], rank=2)
    assert rows[0][3] == "g_rank3_offtarget_mm0"   # rank is 1-based in the name


def test_amplicon_bed_rows_includes_placement_tag():
    rows = pd.amplicon_bed_rows("g", [("c", 0, 100, 100, 0)], rank=0,
                                placement="upstream->internal")
    assert "upstream->internal" in rows[0][3]


def test_amplicon_bed_rows_empty_amplicons():
    assert pd.amplicon_bed_rows("g", []) == []


def test_write_bed_sorts_and_writes_header(tmp_path):
    rows = [
        ["chr2", 10, 50, "a_offtarget_mm0", 1000, "+"],
        ["chr1", 300, 400, "b_offtarget_mm1", 750, "+"],
        ["chr1", 10, 50, "c_ontarget_mm0", 1000, "+"],
    ]
    path = str(tmp_path / "amplicons.bed")
    n = pd.write_bed(path, rows)
    assert n == 3
    with open(path) as fh:
        lines = fh.read().splitlines()
    assert lines[0].startswith("track name=")
    data = lines[1:]
    # Sorted by (chrom, start, end): chr1/10, chr1/300, chr2/10.
    assert [ln.split("\t")[0:2] for ln in data] == [
        ["chr1", "10"], ["chr1", "300"], ["chr2", "10"]]
    # Tab-separated BED6.
    assert data[0].split("\t") == ["chr1", "10", "50", "c_ontarget_mm0", "1000", "+"]


def test_write_bed_empty_still_writes_header(tmp_path):
    path = str(tmp_path / "empty.bed")
    assert pd.write_bed(path, []) == 0
    with open(path) as fh:
        assert fh.read().startswith("track name=")


def test_in_silico_pcr_amplicon_carries_mismatch_count():
    fwd = "AAAGGGCTACGTTGCAATCG"
    rev = "CCCTTTGCAGTACTTAGGAC"
    rev_rc = str(Seq(rev).reverse_complement())
    spacer = "A" * 100
    fwd_5mm = "CC" + fwd[2:]           # off-target: 2 mismatches, 3' seed intact
    template = ("GG" + fwd + spacer + rev_rc
                + ("T" * 300) + fwd_5mm + spacer + rev_rc + "GG")
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=50, max_product=200,
                           seed_len=12, max_mismatches=2)
    mismatches = sorted(a[4] for a in res["amplicons"])
    assert mismatches == [0, 2]      # perfect site + 2-mismatch off-target


def test_in_silico_pcr_unequal_primer_lengths_large_product():
    """Regression for the early-break bug: with primers of different lengths a
    large but in-window amplicon must not be skipped by the sort/break prune."""
    fwd = "AAAGGGCTACGTTGCAATCG"            # 20 bp
    rev = "CCCTTTGCAGTACTTAGGACAGTCA"       # 25 bp
    rev_rc = str(Seq(rev).reverse_complement())
    spacer = "A" * 400
    template = "GG" + fwd + spacer + rev_rc + "GG"
    expected = len(fwd) + len(spacer) + len(rev)
    res = pd.in_silico_pcr(fwd, rev, _genome_from_str(template),
                           min_product=50, max_product=expected)
    assert res["count"] == 1
    assert res["amplicons"][0][3] == expected


# --------------------------------------------------------------------------- #
# Primer design end-to-end on a synthetic template
# --------------------------------------------------------------------------- #
def _random_seq(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_design_primers_success():
    template = _random_seq(800, seed=42)
    params = pd.PrimerParams(product_min=100, product_max=400)
    r = pd.design_primers_for_sequence("synthetic", template, params)
    assert r["status"] == "OK"
    assert r["forward"] and r["reverse"]
    assert r["product_size"] is not None
    assert r["fwd_tm"] is not None
    # reported Tm should be within the requested window (consistency check)
    assert params.min_tm - 5 <= r["fwd_tm"] <= params.max_tm + 5
    assert r["heterodimer_tm"] is not None


def test_design_primer_candidates_ranked_alternates():
    """Requesting N candidates returns up to N distinct, rank-ordered pairs."""
    template = _random_seq(1500, seed=11)
    params = pd.PrimerParams(product_min=100, product_max=600)
    cands = pd.design_primer_candidates("synthetic", template, params, num_candidates=3)
    assert 1 <= len(cands) <= 3
    # Ranks are 0,1,2... in order and all report OK with a primer pair.
    assert [c["rank"] for c in cands] == list(range(len(cands)))
    assert all(c["status"] == "OK" and c["forward"] and c["reverse"] for c in cands)
    if len(cands) > 1:
        # Alternates must be genuinely different oligos, not duplicates.
        pairs = {(c["forward"], c["reverse"]) for c in cands}
        assert len(pairs) == len(cands)


def test_design_primers_for_sequence_is_first_candidate():
    """The single-result wrapper must equal the top-ranked candidate."""
    template = _random_seq(900, seed=12)
    params = pd.PrimerParams(product_min=100, product_max=400)
    single = pd.design_primers_for_sequence("g", template, params)
    first = pd.design_primer_candidates("g", template, params, num_candidates=2)[0]
    assert single["forward"] == first["forward"]
    assert single["reverse"] == first["reverse"]
    assert single["rank"] == 0


def test_design_for_gene_num_candidates_multiplies_rows():
    genome, gene = _placement_genome()
    params = pd.PrimerParams(product_min=100, product_max=500)
    results = pd.design_for_gene(genome, gene, params, flank_size=150,
                                 num_candidates=3)
    assert len(results) >= 1
    assert all(r["placement"] == "internal->internal" for r in results)
    assert [r["rank"] for r in results] == list(range(len(results)))


def test_design_primers_template_too_short():
    r = pd.design_primers_for_sequence("x", "ATGC", pd.PrimerParams())
    assert "too short" in r["status"].lower()


def test_design_primers_shorter_than_product():
    template = _random_seq(80, seed=7)
    params = pd.PrimerParams(product_min=200, product_max=400)
    r = pd.design_primers_for_sequence("x", template, params)
    assert "shorter than min product" in r["status"].lower()


def test_params_validation():
    bad = pd.PrimerParams(min_size=30, opt_size=20, max_size=25)
    assert bad.validate()
    good = pd.PrimerParams()
    assert good.validate() == []


# --------------------------------------------------------------------------- #
# Primer placement (upstream / internal / downstream control)
# --------------------------------------------------------------------------- #
def _placement_genome():
    """Genome with a + strand gene at 200..800 on a 1000 bp contig."""
    seq = _random_seq(1000, seed=99)
    gene = {"chrom": "chr1", "gene_name": "g", "locus_tag": "L1",
            "start": 200, "end": 800, "strand": 1}
    return {"chr1": SeqRecord(Seq(seq), id="chr1")}, gene


def test_placement_combos():
    assert pd.placement_combos("internal") == [("internal", "internal")]
    assert pd.placement_combos("flanking") == [("upstream", "downstream")]
    assert pd.placement_combos(("upstream", "internal")) == [("upstream", "internal")]
    # 'all' = valid permutations where forward is not 3' of reverse
    combos = pd.placement_combos("all")
    assert ("upstream", "downstream") in combos
    assert ("downstream", "upstream") not in combos
    assert len(combos) == 6


def test_build_gene_template_regions():
    genome, gene = _placement_genome()
    template, regions = pd.build_gene_template(genome, gene, flank_size=150)
    assert regions["upstream"] == (0, 150)
    assert regions["internal"] == (150, 600)
    assert regions["downstream"] == (750, 150)
    assert len(template) == 900


def test_design_flanking_places_primers_in_flanks():
    genome, gene = _placement_genome()
    params = pd.PrimerParams(product_min=100, product_max=2000)
    results = pd.design_for_gene(genome, gene, params, flank_size=150, mode="flanking")
    assert len(results) == 1
    r = results[0]
    assert r["placement"] == "upstream->downstream"
    assert r["status"] == "OK"
    # Forward primer should sit in the upstream flank, reverse in downstream flank;
    # i.e. the amplicon must be larger than the gene body (600 bp).
    assert r["product_size"] > 600


def test_design_all_permutations():
    genome, gene = _placement_genome()
    params = pd.PrimerParams(product_min=80, product_max=2000)
    results = pd.design_for_gene(genome, gene, params, flank_size=150, mode="all")
    assert len(results) == 6
    placements = {r["placement"] for r in results}
    assert "upstream->downstream" in placements
    assert "internal->internal" in placements


def test_design_internal_is_default_single_row():
    genome, gene = _placement_genome()
    params = pd.PrimerParams(product_min=100, product_max=500)
    results = pd.design_for_gene(genome, gene, params, flank_size=150)
    assert len(results) == 1
    assert results[0]["placement"] == "internal->internal"


def test_placement_region_too_small():
    genome, gene = _placement_genome()
    params = pd.PrimerParams()
    # No flank -> upstream region is empty, so a flanking design must fail clearly.
    results = pd.design_for_gene(genome, gene, params, flank_size=0, mode="flanking")
    assert "region too small" in results[0]["status"].lower()


def test_p3_template_preserves_length():
    assert len(pd._p3_template("acgtN-x ")) == len("acgtN-x ")
    assert pd._p3_template("acgt") == "ACGT"
    assert "N" in pd._p3_template("ACGT-ACGT")


def test_result_row_matches_columns():
    template = _random_seq(800, seed=3)
    r = pd.design_primers_for_sequence("g", template, pd.PrimerParams())
    r["specificity"] = "Not tested"
    row = pd.result_to_row(r)
    assert len(row) == len(pd.RESULT_COLUMNS)


# --------------------------------------------------------------------------- #
# Ranked alternate candidates (num_return)
# --------------------------------------------------------------------------- #
def test_design_candidates_returns_multiple_ranked():
    """Requesting N candidates returns up to N distinct, rank-ordered pairs."""
    template = _random_seq(800, seed=42)
    params = pd.PrimerParams(product_min=100, product_max=400, num_return=3)
    cands = pd.design_primer_candidates("g", template, params)
    assert 2 <= len(cands) <= 3              # this template yields several
    assert all(c["status"] == "OK" for c in cands)
    # Ranks are 0-based, contiguous and ordered.
    assert [c["rank"] for c in cands] == list(range(len(cands)))
    # The alternates are genuinely different primer pairs, not duplicates.
    pairs = {(c["forward"], c["reverse"]) for c in cands}
    assert len(pairs) == len(cands)
    # Rank in the output row is 1-based for human readers.
    assert pd.result_to_row(cands[0])[2] == 1
    assert pd.result_to_row(cands[1])[2] == 2


def test_design_candidates_default_is_single():
    """Default params (num_return=1) yield exactly one rank-0 candidate."""
    template = _random_seq(800, seed=42)
    cands = pd.design_primer_candidates("g", template,
                                        pd.PrimerParams(product_min=100, product_max=400))
    assert len(cands) == 1 and cands[0]["rank"] == 0


def test_design_primers_for_sequence_returns_best():
    """The single-pair wrapper returns the rank-0 candidate even when
    num_return > 1, so existing callers are unaffected."""
    template = _random_seq(800, seed=42)
    params = pd.PrimerParams(product_min=100, product_max=400, num_return=5)
    best = pd.design_primers_for_sequence("g", template, params)
    cands = pd.design_primer_candidates("g", template, params)
    assert best["rank"] == 0
    assert (best["forward"], best["reverse"]) == (cands[0]["forward"], cands[0]["reverse"])


def test_design_candidates_failure_returns_single_row():
    """A template too short to design still returns exactly one explanatory row."""
    params = pd.PrimerParams(num_return=4)
    cands = pd.design_primer_candidates("x", "ATGC", params)
    assert len(cands) == 1
    assert "too short" in cands[0]["status"].lower()


def test_design_for_gene_emits_candidate_rows():
    """num_return > 1 produces multiple rows per placement for one gene."""
    genome, gene = _placement_genome()
    params = pd.PrimerParams(product_min=100, product_max=500, num_return=3)
    results = pd.design_for_gene(genome, gene, params, flank_size=150)
    ok = [r for r in results if r["status"] == "OK"]
    assert len(ok) >= 2                       # internal placement, several alternates
    assert all(r["placement"] == "internal->internal" for r in ok)
    assert {r["rank"] for r in ok} == set(range(len(ok)))


# --------------------------------------------------------------------------- #
# Quality warnings (Tm balance + secondary structure)
# --------------------------------------------------------------------------- #
def test_design_reports_tm_diff():
    template = _random_seq(800, seed=42)
    r = pd.design_primers_for_sequence("g", template,
                                       pd.PrimerParams(product_min=100, product_max=400))
    assert r["status"] == "OK"
    assert r["tm_diff"] == round(abs(r["fwd_tm"] - r["rev_tm"]), 2)


def test_primer_warnings_flags_problems():
    bad = {
        "fwd_tm": 60.0, "rev_tm": 68.0,          # ΔTm 8°C
        "fwd_hairpin_tm": 55.0, "rev_hairpin_tm": 10.0,
        "fwd_homodimer_tm": 10.0, "rev_homodimer_tm": 10.0,
        "heterodimer_tm": 50.0,
    }
    warns = pd.primer_warnings(bad)
    joined = " ".join(warns)
    assert any("ΔTm" in w for w in warns)
    assert "hairpin" in joined and "hetero-dimer" in joined


def test_primer_warnings_clean_pair_is_empty():
    good = {
        "fwd_tm": 60.0, "rev_tm": 60.5,
        "fwd_hairpin_tm": 20.0, "rev_hairpin_tm": 18.0,
        "fwd_homodimer_tm": 5.0, "rev_homodimer_tm": 5.0,
        "heterodimer_tm": 12.0,
    }
    assert pd.primer_warnings(good) == []


def test_primer_warnings_tolerates_missing_values():
    # None values (failed Tm/structure calc) must not raise.
    assert pd.primer_warnings({}) == []
    assert pd.primer_warnings({"fwd_tm": None, "rev_tm": 60.0}) == []


# --------------------------------------------------------------------------- #
# Composite quality score + quality re-ranking
# --------------------------------------------------------------------------- #
def test_quality_score_clean_pair_is_high():
    """A near-ideal pair (Tm == opt, balanced, 50% GC, no structure) scores ~100."""
    clean = {
        "forward": "ACGT", "reverse": "ACGT",
        "fwd_tm": 60.0, "rev_tm": 60.0, "tm_diff": 0.0,
        "fwd_gc": 50.0, "rev_gc": 50.0,
        "fwd_hairpin_tm": 10.0, "rev_hairpin_tm": 10.0,
        "fwd_homodimer_tm": 10.0, "rev_homodimer_tm": 10.0,
        "heterodimer_tm": 10.0,
    }
    assert pd.quality_score(clean, pd.PrimerParams(opt_tm=60.0)) == 100.0


def test_quality_score_penalises_problems():
    """A pair with Tm imbalance and a stable hetero-dimer scores lower."""
    good = {
        "forward": "ACGT", "reverse": "ACGT",
        "fwd_tm": 60.0, "rev_tm": 60.0, "tm_diff": 0.0,
        "fwd_gc": 50.0, "rev_gc": 50.0,
        "fwd_hairpin_tm": 10.0, "rev_hairpin_tm": 10.0,
        "fwd_homodimer_tm": 10.0, "rev_homodimer_tm": 10.0,
        "heterodimer_tm": 10.0,
    }
    bad = dict(good, fwd_tm=66.0, tm_diff=6.0, heterodimer_tm=55.0)
    g = pd.quality_score(good, pd.PrimerParams(opt_tm=60.0))
    b = pd.quality_score(bad, pd.PrimerParams(opt_tm=60.0))
    assert b < g
    assert 0.0 <= b <= 100.0


def test_quality_score_none_for_failure_row():
    """A result with no primer pair has no score (blank column)."""
    assert pd.quality_score({"forward": None, "reverse": None}) is None
    assert pd.quality_score(pd._new_result("g", "internal")) is None


def test_quality_score_clamped_to_zero():
    """A catastrophic pair never scores below zero."""
    awful = {
        "forward": "ACGT", "reverse": "ACGT",
        "fwd_tm": 90.0, "rev_tm": 30.0, "tm_diff": 60.0,
        "fwd_gc": 100.0, "rev_gc": 0.0,
        "fwd_hairpin_tm": 80.0, "rev_hairpin_tm": 80.0,
        "fwd_homodimer_tm": 80.0, "rev_homodimer_tm": 80.0,
        "heterodimer_tm": 80.0,
    }
    assert pd.quality_score(awful, pd.PrimerParams(opt_tm=60.0)) == 0.0


def test_designed_candidate_carries_quality_score():
    """End-to-end: a designed OK pair has a numeric quality score in range."""
    template = _random_seq(800, seed=42)
    r = pd.design_primers_for_sequence("g", template,
                                       pd.PrimerParams(product_min=100, product_max=400))
    assert r["status"] == "OK"
    assert isinstance(r["quality_score"], float)
    assert 0.0 <= r["quality_score"] <= 100.0
    # The score must appear in the flattened CSV row.
    r["specificity"] = "Not tested"
    row = pd.result_to_row(r)
    assert len(row) == len(pd.RESULT_COLUMNS)
    assert row[pd.RESULT_COLUMNS.index("Quality Score")] == r["quality_score"]


def test_rank_by_quality_orders_by_score_desc():
    """rank_by='quality' returns candidates sorted by descending quality score
    with contiguous 0-based ranks; the same pairs, just reordered."""
    template = _random_seq(1500, seed=11)
    params = pd.PrimerParams(product_min=100, product_max=600, num_return=4)
    default = pd.design_primer_candidates("g", template, params)
    byq = pd.design_primer_candidates("g", template, params, rank_by="quality")
    if len(byq) > 1:
        scores = [c["quality_score"] for c in byq]
        assert scores == sorted(scores, reverse=True)
    assert [c["rank"] for c in byq] == list(range(len(byq)))
    # Re-ranking is a permutation: identical set of primer pairs either way.
    assert ({(c["forward"], c["reverse"]) for c in default}
            == {(c["forward"], c["reverse"]) for c in byq})


def test_rank_by_quality_can_promote_cleaner_pair():
    """The rank-0 pair under quality ranking has the best score of the set."""
    template = _random_seq(1500, seed=11)
    params = pd.PrimerParams(product_min=100, product_max=600, num_return=5)
    byq = pd.design_primer_candidates("g", template, params, rank_by="quality")
    best = max(c["quality_score"] for c in byq)
    assert byq[0]["quality_score"] == best


def test_rank_by_invalid_raises():
    template = _random_seq(800, seed=42)
    with pytest.raises(ValueError):
        pd.design_primer_candidates("g", template,
                                    pd.PrimerParams(product_min=100, product_max=400),
                                    rank_by="bogus")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
