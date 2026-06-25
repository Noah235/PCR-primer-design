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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
