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

    return legacy["elapsed"], fixed["elapsed"], len(primer_pairs)


def bench_design(templates):
    params = pd.PrimerParams(product_min=100, product_max=600)
    ok = 0
    with timer() as t:
        for i, tpl in enumerate(templates):
            r = pd.design_primers_for_sequence(f"g{i}", tpl, params)
            ok += r["status"] == "OK"
    return t["elapsed"], ok, len(templates)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-bp", type=int, default=2_000_000)
    ap.add_argument("--pairs", type=int, default=20, help="primer pairs for specificity bench")
    ap.add_argument("--genes", type=int, default=30, help="templates for design bench")
    ap.add_argument("--template-bp", type=int, default=1200)
    ap.add_argument("-o", "--output", default="benchmark_results.md")
    args = ap.parse_args()

    print(f"Building synthetic genome: {args.genome_bp:,} bp ...")
    genome_seq, genome = make_genome(args.genome_bp)
    fasta_path = "_bench_genome.fasta"
    SeqIO.write(genome["chr1"], fasta_path, "fasta")

    print(f"Specificity benchmark: {args.pairs} primer pairs ...")
    legacy_t, fixed_t, n_pairs = bench_specificity(genome_seq, genome, fasta_path, args.pairs)

    print(f"Design benchmark: {args.genes} templates ...")
    templates = make_templates(genome_seq, args.genes, args.template_bp)
    design_t, ok, n_t = bench_design(templates)

    speedup = legacy_t / fixed_t if fixed_t else float("inf")

    lines = [
        "# Benchmark Results",
        "",
        f"- Genome size: **{args.genome_bp:,} bp**",
        f"- Python/primer3 in-process; single thread.",
        "",
        "## Specificity search (in-silico PCR)",
        "",
        "| Implementation | Pairs | Total time | Per pair |",
        "| --- | ---: | ---: | ---: |",
        f"| Legacy (reload genome per pair, 1 orientation) | {n_pairs} | "
        f"{legacy_t:.3f} s | {legacy_t / n_pairs * 1000:.1f} ms |",
        f"| Fixed (preloaded genome, 2 orientations) | {n_pairs} | "
        f"{fixed_t:.3f} s | {fixed_t / n_pairs * 1000:.1f} ms |",
        "",
        f"**Speedup: {speedup:.1f}x** (the fixed version also checks both primer "
        "orientations, which the legacy version missed).",
        "",
        "## Primer design throughput",
        "",
        "| Templates | Successful | Total time | Per template |",
        "| ---: | ---: | ---: | ---: |",
        f"| {n_t} | {ok} | {design_t:.3f} s | {design_t / n_t * 1000:.1f} ms |",
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
