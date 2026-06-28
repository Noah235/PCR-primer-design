#!/usr/bin/env python3
"""Headless command-line interface for batch PCR primer design.

Examples
--------
Genome + GFF3, design primers for three genes with specificity testing::

    python primer_cli.py genome \\
        --genome genome.fasta --gff annotation.gff3 \\
        --genes sulA,opgH,galU --specificity -o primers.csv

CDS-only mode for every record in a CDS FASTA::

    python primer_cli.py cds --cds cds.fasta -o primers.csv

The CLI shares 100%% of its logic with the GUI (both call :mod:`primer_design`),
so results are identical between the two front-ends.
"""

import argparse
import csv
import logging
import sys

import primer_design as pd


def _add_param_args(p):
    g = p.add_argument_group("primer parameters")
    d = pd.PrimerParams()
    g.add_argument("--min-size", type=int, default=d.min_size)
    g.add_argument("--opt-size", type=int, default=d.opt_size)
    g.add_argument("--max-size", type=int, default=d.max_size)
    g.add_argument("--min-tm", type=float, default=d.min_tm)
    g.add_argument("--opt-tm", type=float, default=d.opt_tm)
    g.add_argument("--max-tm", type=float, default=d.max_tm)
    g.add_argument("--min-gc", type=float, default=d.min_gc)
    g.add_argument("--max-gc", type=float, default=d.max_gc)
    g.add_argument("--product-min", type=int, default=d.product_min)
    g.add_argument("--product-max", type=int, default=d.product_max)
    g.add_argument("--gc-clamp", type=int, default=d.gc_clamp)
    g.add_argument("--num-return", type=int, default=d.num_return,
                   help="ranked primer-pair candidates to report per template "
                        "(1 = best only; each extra pair is an additional row)")


def _params_from_args(a) -> pd.PrimerParams:
    params = pd.PrimerParams(
        min_size=a.min_size, opt_size=a.opt_size, max_size=a.max_size,
        min_tm=a.min_tm, opt_tm=a.opt_tm, max_tm=a.max_tm,
        min_gc=a.min_gc, max_gc=a.max_gc,
        product_min=a.product_min, product_max=a.product_max,
        gc_clamp=a.gc_clamp, num_return=a.num_return,
    )
    problems = params.validate()
    if problems:
        sys.exit("Invalid parameters: " + "; ".join(problems))
    return params


def _write(output_csv, banner, rows):
    with open(output_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([banner])
        w.writerow(pd.RESULT_COLUMNS)
        for row in rows:
            w.writerow(row)


def _resolve_placement(a):
    """Turn CLI --placement / --fwd-region / --rev-region into a placement mode."""
    if a.placement == "custom":
        if not a.fwd_region or not a.rev_region:
            sys.exit("--placement custom requires --fwd-region and --rev-region "
                     f"(choose from {', '.join(pd.PLACEMENT_REGIONS)}).")
        if a.fwd_region not in pd.PLACEMENT_REGIONS or a.rev_region not in pd.PLACEMENT_REGIONS:
            sys.exit(f"regions must be one of: {', '.join(pd.PLACEMENT_REGIONS)}")
        return (a.fwd_region, a.rev_region)
    return a.placement


def run_genome(a):
    params = _params_from_args(a)
    placement = _resolve_placement(a)
    genome = pd.load_genome(a.genome)
    all_genes = pd.parse_gff3_full(a.gff)
    gene_list, _found, not_found = pd.filter_genes_by_names(all_genes, a.genes or "")
    if not_found:
        logging.warning("Genes not found: %s", ", ".join(not_found))
    if not gene_list:
        sys.exit("No genes to process.")
    logging.info("Designing primers for %d gene(s), placement=%s ...",
                 len(gene_list), placement)

    prepared_genome = pd.prepare_genome(genome) if a.specificity else None
    rows, n_ok = [], 0
    for gene in gene_list:
        for r in pd.design_for_gene(genome, gene, params, a.flank, mode=placement):
            if a.specificity and r["status"] == "OK":
                max_prod = max(params.product_max, (r["product_size"] or 0) + 500)
                spec = pd.in_silico_pcr(r["forward"], r["reverse"], prepared_genome,
                                        50, max_prod,
                                        seed_len=a.seed_len,
                                        max_mismatches=a.max_mismatches)
                r["specificity"] = pd.specificity_label(spec)
            else:
                r["specificity"] = "Not tested"
            n_ok += r["status"] == "OK"
            rows.append(pd.result_to_row(r))

    if a.specificity:
        spec_desc = (f"specificity=seed{a.seed_len}/mm{a.max_mismatches}"
                     if a.seed_len else "specificity=exact")
    else:
        spec_desc = "specificity=off"
    _write(a.output, pd.params_summary(
        params,
        f"placement={placement}, flank={a.flank}, genes={len(gene_list)}, {spec_desc}"), rows)
    logging.info("%d primer pair(s) designed across %d genes -> %s",
                 n_ok, len(gene_list), a.output)


def run_cds(a):
    params = _params_from_args(a)
    all_cds = pd.load_cds_sequences(a.cds)
    names = pd.parse_target_names(a.genes or "")
    if names:
        wanted = set(names)
        cds = {k: v for k, v in all_cds.items() if k in wanted}
        for n in (n for n in names if n not in all_cds):
            logging.warning("CDS not found: %s", n)
    else:
        cds = all_cds
    if not cds:
        sys.exit("No CDS to process.")
    logging.info("Designing primers for %d CDS...", len(cds))

    rows, n_ok = [], 0
    for name, info in cds.items():
        for r in pd.design_primer_candidates(name, info["sequence"], params):
            r["specificity"] = "N/A (CDS mode)"
            n_ok += r["status"] == "OK"
            rows.append(pd.result_to_row(r))

    _write(a.output, pd.params_summary(params, f"mode=CDS, cds={len(cds)}"), rows)
    logging.info("%d/%d CDS had suitable primers -> %s", n_ok, len(cds), a.output)


def build_parser():
    p = argparse.ArgumentParser(description="Batch PCR primer design (headless).")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = p.add_subparsers(dest="command", required=True)

    g = sub.add_parser("genome", help="genome FASTA + GFF3 pipeline")
    g.add_argument("--genome", required=True)
    g.add_argument("--gff", required=True)
    g.add_argument("--genes", default="", help="comma/line-separated names; empty = all")
    g.add_argument("--specificity", action="store_true", help="run in-silico PCR check")
    g.add_argument("--seed-len", type=int, default=0,
                   help="3'-anchor seed length for mismatch-tolerant specificity "
                        "(0 = exact match; ~12 is a sensible value)")
    g.add_argument("--max-mismatches", type=int, default=0,
                   help="max mismatches in the primer 5' of the seed "
                        "(only used when --seed-len > 0)")
    g.add_argument("--flank", type=int, default=200,
                   help="flank size (bp) extracted up/downstream of each gene")
    g.add_argument("--placement", default="internal",
                   choices=["internal", "flanking", "all", "custom"],
                   help="where to place primers relative to the gene "
                        "(custom: use --fwd-region/--rev-region)")
    g.add_argument("--fwd-region", choices=list(pd.PLACEMENT_REGIONS),
                   help="forward-primer region for --placement custom")
    g.add_argument("--rev-region", choices=list(pd.PLACEMENT_REGIONS),
                   help="reverse-primer region for --placement custom")
    g.add_argument("-o", "--output", default="primers.csv")
    _add_param_args(g)
    g.set_defaults(func=run_genome)

    c = sub.add_parser("cds", help="CDS-FASTA-only pipeline")
    c.add_argument("--cds", required=True)
    c.add_argument("--genes", default="", help="comma/line-separated names; empty = all")
    c.add_argument("-o", "--output", default="primers.csv")
    _add_param_args(c)
    c.set_defaults(func=run_cds)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )
    args.func(args)


if __name__ == "__main__":
    main()
