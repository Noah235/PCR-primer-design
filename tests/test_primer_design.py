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
    chrom, start, end, size = res["amplicons"][0]
    assert size == len(fwd) + len(spacer) + len(rev)


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


def test_specificity_label():
    assert pd.specificity_label({"count": 0}) == "No amplicons found"
    assert pd.specificity_label({"count": 1}).startswith("Specific")
    assert pd.specificity_label({"count": 3}).startswith("Non-specific")
    assert pd.specificity_label({"count": -1, "error": "bad"}) == "bad"


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


def test_result_row_matches_columns():
    template = _random_seq(800, seed=3)
    r = pd.design_primers_for_sequence("g", template, pd.PrimerParams())
    r["specificity"] = "Not tested"
    row = pd.result_to_row(r)
    assert len(row) == len(pd.RESULT_COLUMNS)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
