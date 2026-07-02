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
                   help="ranked primer pairs to report per target (1 = best only)")
    g.add_argument("--rank-by-quality", action="store_true",
                   help="re-order candidates by the composite quality score "
                        "(accounts for hairpin/dimer Tm) instead of Primer3's "
                        "default ranking, so a cleaner pair can become rank 1")


def _params_from_args(a) -> pd.PrimerParams:
    params = pd.PrimerParams(
        min_size=a.min_size, opt_size=a.opt_size, max_size=a.max_size,
        min_tm=a.min_tm, opt_tm=a.opt_tm, max_tm=a.max_tm,
        min_gc=a.min_gc, max_gc=a.max_gc,
        product_min=a.product_min, product_max=a.product_max,
        gc_clamp=a.gc_clamp, num_return=max(1, a.num_return),
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
    if a.bed and not a.specificity:
        sys.exit("--bed requires --specificity (amplicon locations come from the "
                 "in-silico PCR specificity search).")
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
    rank_by = "quality" if a.rank_by_quality else "primer3"
    rows, n_ok, bed_rows = [], 0, []
    for gene in gene_list:
        target = (gene["chrom"], gene["start"], gene["end"])
        for r in pd.design_for_gene(genome, gene, params, a.flank, mode=placement,
                                    num_candidates=params.num_return, rank_by=rank_by):
            if a.specificity and r["status"] == "OK":
                max_prod = max(params.product_max, (r["product_size"] or 0) + 500)
                spec = pd.in_silico_pcr(r["forward"], r["reverse"], prepared_genome,
                                        50, max_prod,
                                        seed_len=a.seed_len,
                                        max_mismatches=a.max_mismatches)
                r["specificity"] = pd.specificity_label(spec)
                if a.bed:
                    bed_rows.extend(pd.amplicon_bed_rows(
                        gene["gene_name"] or gene["locus_tag"], spec["amplicons"],
                        rank=r.get("rank", 0), target=target,
                        placement=r.get("placement", "")))
            else:
                r["specificity"] = "Not tested"
            n_ok += r["status"] == "OK"
            rows.append(pd.result_to_row(r))

    if a.bed:
        n_bed = pd.write_bed(a.bed, bed_rows)
        logging.info("%d predicted amplicon location(s) -> %s", n_bed, a.bed)

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

    rank_by = "quality" if a.rank_by_quality else "primer3"
    rows, n_ok = [], 0
    for name, info in cds.items():
        for r in pd.design_primer_candidates(name, info["sequence"], params,
                                             num_candidates=params.num_return,
                                             rank_by=rank_by):
            r["specificity"] = "N/A (CDS mode)"
            n_ok += r["status"] == "OK"
            rows.append(pd.result_to_row(r))

    _write(a.output, pd.params_summary(params, f"mode=CDS, cds={len(cds)}"), rows)
    logging.info("%d primer pair(s) designed across %d CDS -> %s",
                 n_ok, len(cds), a.output)


def run_check(a):
    """In-silico QC of user-supplied primers (no template design).

    Reports each pair's Tm, GC%, secondary-structure Tm, composite quality score
    and warnings, and — when a genome is given — runs the accurate
    two-orientation (optionally mismatch-tolerant) in-silico PCR specificity
    check and reports the predicted product size and off-target count.
    """
    params = _params_from_args(a)
    pairs = []
    if a.primers:
        pairs.extend(pd.read_primer_pairs(a.primers))
    if a.forward and a.reverse:
        pairs.append((a.name or f"pair{len(pairs) + 1}", a.forward, a.reverse))
    if not pairs:
        sys.exit("Provide primers via --forward/--reverse or --primers FILE.")
    if a.bed and not a.genome:
        sys.exit("--bed requires --genome (amplicon locations come from the "
                 "in-silico PCR search against a genome).")

    prepared = None
    if a.genome:
        prepared = pd.prepare_genome(pd.load_genome(a.genome))
        logging.info("Genome loaded; checking %d primer pair(s)...", len(pairs))
    else:
        logging.info("No genome given: reporting primer QC only for %d pair(s).",
                     len(pairs))

    # Generous search window so a pair's own (possibly large) product is found
    # while still catching genuine off-target products a few kb apart.
    min_prod, max_prod = 50, max(params.product_max, 5000)
    rows, bed_rows, n_specific = [], [], 0
    for name, fwd, rev in pairs:
        r = pd.analyze_primer_pair(name, fwd, rev, params)
        if prepared is not None and r["status"] == "OK":
            spec = pd.in_silico_pcr(r["forward"], r["reverse"], prepared,
                                    min_prod, max_prod,
                                    seed_len=a.seed_len,
                                    max_mismatches=a.max_mismatches)
            r["specificity"] = pd.specificity_label(spec)
            # A single predicted amplicon is unambiguously the product size.
            if spec.get("count") == 1:
                r["product_size"] = spec["amplicons"][0][3]
            n_specific += spec.get("count") == 1
            if a.bed:
                bed_rows.extend(pd.amplicon_bed_rows(
                    name, spec["amplicons"], rank=0, target=None,
                    placement="user"))
        elif prepared is None:
            r["specificity"] = "Not tested (no genome)"
        rows.append(pd.result_to_row(r))

    if a.bed:
        n_bed = pd.write_bed(a.bed, bed_rows)
        logging.info("%d predicted amplicon location(s) -> %s", n_bed, a.bed)

    spec_desc = "no genome"
    if prepared is not None:
        spec_desc = (f"specificity=seed{a.seed_len}/mm{a.max_mismatches}"
                     if a.seed_len else "specificity=exact")
    _write(a.output, pd.params_summary(
        params, f"mode=check, pairs={len(pairs)}, {spec_desc}"), rows)
    if prepared is not None:
        logging.info("%d/%d pair(s) specific (single predicted amplicon) -> %s",
                     n_specific, len(pairs), a.output)
    else:
        logging.info("%d primer pair(s) QC'd -> %s", len(pairs), a.output)


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
    g.add_argument("--bed", default=None,
                   help="also write predicted amplicon locations (intended + "
                        "off-target) to this BED file for a genome browser "
                        "(requires --specificity)")
    _add_param_args(g)
    g.set_defaults(func=run_genome)

    c = sub.add_parser("cds", help="CDS-FASTA-only pipeline")
    c.add_argument("--cds", required=True)
    c.add_argument("--genes", default="", help="comma/line-separated names; empty = all")
    c.add_argument("-o", "--output", default="primers.csv")
    _add_param_args(c)
    c.set_defaults(func=run_cds)

    ch = sub.add_parser(
        "check",
        help="QC / in-silico-PCR check of primers you already have")
    ch.add_argument("--forward", help="forward primer sequence (5'->3')")
    ch.add_argument("--reverse", help="reverse primer sequence (5'->3')")
    ch.add_argument("--name", default=None, help="label for a single pair")
    ch.add_argument("--primers", default=None,
                    help="batch file of 'name forward reverse' lines "
                         "(comma/tab/space separated; '#' comments ignored)")
    ch.add_argument("--genome", default=None,
                    help="optional genome FASTA to run the specificity check "
                         "against (omit for primer QC only)")
    ch.add_argument("--seed-len", type=int, default=0,
                    help="3'-anchor seed length for mismatch-tolerant "
                         "specificity (0 = exact; ~12 is sensible)")
    ch.add_argument("--max-mismatches", type=int, default=0,
                    help="max 5' mismatches (only used when --seed-len > 0)")
    ch.add_argument("--bed", default=None,
                    help="write predicted amplicon locations to BED "
                         "(requires --genome)")
    ch.add_argument("-o", "--output", default="primer_check.csv")
    _add_param_args(ch)
    ch.set_defaults(func=run_check)
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
