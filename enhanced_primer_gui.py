"""Tkinter GUI for the PCR primer-design suite.

All sequence/primer logic lives in :mod:`primer_design`; this file is a thin
presentation layer. Importing this module no longer launches a window — the GUI
only starts when run as ``python enhanced_primer_gui.py`` (or via
:func:`main`), so the logic can be imported and tested headlessly.
"""

import csv
import os

import primer_design as pd


def _run_genome_gff(params, output_csv, genome_path, gff_path, target_names,
                    flank_size, do_specificity, log, placement="internal",
                    seed_len=0, max_mismatches=0):
    """Genome FASTA + GFF3 pipeline. ``log`` is a callable(str)."""
    genome = pd.load_genome(genome_path)
    log("Genome loaded")

    all_genes = pd.parse_gff3_full(gff_path)
    gene_list, found, not_found = pd.filter_genes_by_names(all_genes, target_names)
    if pd.parse_target_names(target_names):
        log(f"Filtered {len(gene_list)} of {len(all_genes)} genes")
        if not_found:
            log(f"Not found: {', '.join(not_found[:5])}")
    else:
        log(f"Processing all {len(gene_list)} genes")

    if not gene_list:
        raise ValueError("No genes found to process!")

    extracted_fasta = os.path.splitext(output_csv)[0] + "_extracted_sequences.fasta"
    pd.extract_sequences_to_fasta(genome, gene_list, extracted_fasta, flank_size=flank_size)
    log("Extracting sequences...")

    # Pre-clean the genome once for specificity testing (big speed win).
    prepared_genome = pd.prepare_genome(genome) if do_specificity else None

    n_ok = 0
    with open(output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([pd.params_summary(
            params, f"flank={flank_size}, placement={placement}, genes={len(gene_list)}")])
        writer.writerow(pd.RESULT_COLUMNS)
        for i, gene in enumerate(gene_list):
            if i % 10 == 0:
                log(f"Processing gene {i + 1}/{len(gene_list)}...")
            for r in pd.design_for_gene(genome, gene, params, flank_size, mode=placement,
                                        num_candidates=params.num_return):
                if do_specificity and r["status"] == "OK":
                    # Window generous enough to include this pair's own amplicon.
                    max_prod = max(params.product_max, (r["product_size"] or 0) + 500)
                    spec = pd.in_silico_pcr(
                        r["forward"], r["reverse"], prepared_genome,
                        min_product=50, max_product=max_prod,
                        seed_len=seed_len, max_mismatches=max_mismatches)
                    r["specificity"] = pd.specificity_label(spec)
                else:
                    r["specificity"] = "Not tested"
                if r["status"] == "OK":
                    n_ok += 1
                writer.writerow(pd.result_to_row(r))

    log(f"Done. {n_ok} primer pair(s) designed across {len(gene_list)} genes")
    log(f"Output: {os.path.abspath(output_csv)}")
    log(f"Sequences: {os.path.abspath(extracted_fasta)}")
    return n_ok, len(gene_list)


def _run_cds(params, output_csv, cds_path, target_names, log):
    """CDS-FASTA-only pipeline."""
    all_cds = pd.load_cds_sequences(cds_path)
    log(f"Loaded {len(all_cds)} CDS sequences")

    names = pd.parse_target_names(target_names)
    if names:
        wanted = set(names)
        filtered = {k: v for k, v in all_cds.items() if k in wanted}
        not_found = [n for n in names if n not in all_cds]
        log(f"Filtered {len(filtered)} of {len(all_cds)} CDS")
        if not_found:
            log(f"Not found: {', '.join(not_found[:5])}")
    else:
        filtered = all_cds
        log(f"Processing all {len(filtered)} CDS")

    if not filtered:
        raise ValueError("No CDS sequences found to process!")

    n_ok = 0
    with open(output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([pd.params_summary(params, f"mode=CDS, cds={len(filtered)}")])
        writer.writerow(pd.RESULT_COLUMNS)
        for i, (name, info) in enumerate(filtered.items()):
            if i % 10 == 0:
                log(f"Processing CDS {i + 1}/{len(filtered)}...")
            for r in pd.design_primer_candidates(name, info["sequence"], params,
                                                 num_candidates=params.num_return):
                r["specificity"] = "N/A (CDS mode)"
                if r["status"] == "OK":
                    n_ok += 1
                writer.writerow(pd.result_to_row(r))

    log(f"Done. {n_ok} primer pair(s) designed across {len(filtered)} CDS")
    log(f"Output: {os.path.abspath(output_csv)}")
    return n_ok, len(filtered)


def main():  # pragma: no cover - interactive GUI
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    root = tk.Tk()
    root.title("🧬 PCR Primer Design Suite")
    root.geometry("960x840")

    # ----- Input mode -----
    mode_frame = ttk.LabelFrame(root, text="Input Mode", padding="10")
    mode_frame.grid(row=0, column=0, columnspan=6, sticky="ew", padx=10, pady=5)
    mode_var = tk.StringVar(value="genome_gff")

    # ----- Files -----
    file_frame = ttk.LabelFrame(root, text="Input Files", padding="10")
    file_frame.grid(row=1, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

    def browse_open(entry):
        path = filedialog.askopenfilename()
        if path:
            entry.delete(0, tk.END)
            entry.insert(0, path)

    genome_label = tk.Label(file_frame, text="Genome FASTA")
    genome_label.grid(row=0, column=0, sticky="w")
    genome_entry = tk.Entry(file_frame, width=55)
    genome_entry.grid(row=0, column=1)
    genome_browse = tk.Button(file_frame, text="Browse", command=lambda: browse_open(genome_entry))
    genome_browse.grid(row=0, column=2)

    gff_label = tk.Label(file_frame, text="GFF3 File")
    gff_label.grid(row=1, column=0, sticky="w")
    gff_entry = tk.Entry(file_frame, width=55)
    gff_entry.grid(row=1, column=1)
    gff_browse = tk.Button(file_frame, text="Browse", command=lambda: browse_open(gff_entry))
    gff_browse.grid(row=1, column=2)

    cds_label = tk.Label(file_frame, text="CDS FASTA")
    cds_entry = tk.Entry(file_frame, width=55)
    cds_browse = tk.Button(file_frame, text="Browse", command=lambda: browse_open(cds_entry))

    tk.Label(file_frame, text="Output CSV").grid(row=2, column=0, sticky="w")
    output_entry = tk.Entry(file_frame, width=55)
    output_entry.insert(0, "primers.csv")
    output_entry.grid(row=2, column=1)
    tk.Button(
        file_frame, text="Browse",
        command=lambda: (lambda p: (output_entry.delete(0, tk.END), output_entry.insert(0, p)) if p else None)(
            filedialog.asksaveasfilename(defaultextension=".csv"))
    ).grid(row=2, column=2)

    # ----- Gene filter -----
    filter_frame = ttk.LabelFrame(root, text="Gene/CDS Selection (case-insensitive)", padding="10")
    filter_frame.grid(row=2, column=0, columnspan=6, sticky="ew", padx=10, pady=5)
    tk.Label(filter_frame, text="Gene names (comma/line separated). '#' lines are ignored.").grid(
        row=0, column=0, sticky="w")
    tk.Label(filter_frame, text="Leave empty to process ALL", fg="blue").grid(row=0, column=1, sticky="w")
    gene_filter_text = scrolledtext.ScrolledText(filter_frame, width=60, height=4)
    gene_filter_text.grid(row=1, column=0, columnspan=2, sticky="ew", pady=5)
    gene_filter_text.insert(tk.END, "# Examples (these comment lines are ignored):\n# sulA, opgH, galU")

    # ----- Parameters -----
    param_frame = ttk.LabelFrame(root, text="Primer Parameters", padding="10")
    param_frame.grid(row=3, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

    def add_entry(label, r, c, default):
        tk.Label(param_frame, text=label).grid(row=r, column=c)
        e = tk.Entry(param_frame, width=6)
        e.insert(0, str(default))
        e.grid(row=r, column=c + 1)
        return e

    defaults = pd.PrimerParams()
    min_size_e = add_entry("Min Size", 0, 0, defaults.min_size)
    opt_size_e = add_entry("Opt Size", 0, 2, defaults.opt_size)
    max_size_e = add_entry("Max Size", 0, 4, defaults.max_size)
    min_tm_e = add_entry("Min Tm", 1, 0, defaults.min_tm)
    opt_tm_e = add_entry("Opt Tm", 1, 2, defaults.opt_tm)
    max_tm_e = add_entry("Max Tm", 1, 4, defaults.max_tm)
    min_gc_e = add_entry("Min GC%", 2, 0, defaults.min_gc)
    max_gc_e = add_entry("Max GC%", 2, 2, defaults.max_gc)
    min_prod_e = add_entry("Min Product", 2, 4, defaults.product_min)
    max_prod_e = add_entry("Max Product", 3, 0, defaults.product_max)
    gc_clamp_e = add_entry("GC Clamp", 3, 2, defaults.gc_clamp)
    num_return_e = add_entry("Candidates", 6, 0, defaults.num_return)
    flank_label = tk.Label(param_frame, text="Flank (bp)")
    flank_label.grid(row=3, column=4)
    flank_e = tk.Entry(param_frame, width=6)
    flank_e.insert(0, "200")
    flank_e.grid(row=3, column=5)

    check_specificity = tk.BooleanVar(value=True)
    tk.Checkbutton(param_frame, text="Test primer specificity",
                   variable=check_specificity).grid(row=4, column=0, columnspan=2, sticky="w")

    # Mismatch-tolerant specificity controls (seed_len=0 -> exact match).
    seed_len_e = add_entry("3' seed (0=exact)", 5, 0, 0)
    max_mm_e = add_entry("Max mismatches", 5, 2, 0)
    # Ranked alternates: how many primer pairs to report per target (1 = best).
    num_return_e = add_entry("Primers/target", 5, 4, defaults.num_return)

    # Primer placement relative to the gene (genome+GFF mode only).
    placement_label = tk.Label(param_frame, text="Primer placement")
    placement_label.grid(row=4, column=2, sticky="e")
    placement_choices = {
        "Internal (both inside gene)": "internal",
        "Flanking (fwd upstream / rev downstream)": "flanking",
        "Fwd upstream / rev internal": ("upstream", "internal"),
        "Fwd internal / rev downstream": ("internal", "downstream"),
        "All permutations": "all",
    }
    placement_var = tk.StringVar(value="Internal (both inside gene)")
    placement_combo = ttk.Combobox(param_frame, textvariable=placement_var, state="readonly",
                                   values=list(placement_choices), width=36)
    placement_combo.grid(row=4, column=3, columnspan=3, sticky="w")

    # ----- Results -----
    results_frame = ttk.LabelFrame(root, text="Results", padding="10")
    results_frame.grid(row=5, column=0, columnspan=6, sticky="ew", padx=10, pady=5)
    results_text = scrolledtext.ScrolledText(results_frame, width=110, height=16)
    results_text.pack(fill=tk.BOTH, expand=True)

    def log(msg):
        results_text.insert(tk.END, msg + "\n")
        results_text.see(tk.END)
        root.update()

    def toggle_mode():
        genome_gff = mode_var.get() == "genome_gff"
        for w in (genome_label, genome_entry, genome_browse, gff_label, gff_entry,
                  gff_browse, flank_label, flank_e, placement_label, placement_combo):
            (w.grid if genome_gff else w.grid_remove)()
        for w in (cds_label, cds_entry, cds_browse):
            (w.grid_remove if genome_gff else w.grid)()
        if not genome_gff:
            cds_label.grid(row=0, column=0, sticky="w")
            cds_entry.grid(row=0, column=1)
            cds_browse.grid(row=0, column=2)

    def build_params():
        params = pd.PrimerParams(
            min_size=int(min_size_e.get()), opt_size=int(opt_size_e.get()),
            max_size=int(max_size_e.get()),
            min_tm=float(min_tm_e.get()), opt_tm=float(opt_tm_e.get()),
            max_tm=float(max_tm_e.get()),
            min_gc=float(min_gc_e.get()), max_gc=float(max_gc_e.get()),
            product_min=int(min_prod_e.get()), product_max=int(max_prod_e.get()),
            gc_clamp=int(gc_clamp_e.get()),
            num_return=max(1, int(num_return_e.get() or 1)),
        )
        problems = params.validate()
        if problems:
            raise ValueError("; ".join(problems))
        return params

    def run_pipeline():
        output_csv = output_entry.get()
        if not output_csv:
            messagebox.showerror("Error", "Please specify an output CSV file.")
            return
        try:
            params = build_params()
        except ValueError as exc:
            messagebox.showerror("Invalid parameters", str(exc))
            return

        results_text.delete("1.0", tk.END)
        log("Starting primer design...")
        target_names = gene_filter_text.get("1.0", tk.END)
        try:
            if mode_var.get() == "genome_gff":
                if not genome_entry.get() or not gff_entry.get():
                    messagebox.showerror("Error", "Select genome FASTA and GFF3 files.")
                    return
                _run_genome_gff(params, output_csv, genome_entry.get(), gff_entry.get(),
                                target_names, int(flank_e.get()), check_specificity.get(), log,
                                placement=placement_choices[placement_var.get()],
                                seed_len=int(seed_len_e.get() or 0),
                                max_mismatches=int(max_mm_e.get() or 0))
            else:
                if not cds_entry.get():
                    messagebox.showerror("Error", "Select a CDS FASTA file.")
                    return
                _run_cds(params, output_csv, cds_entry.get(), target_names, log)
        except Exception as exc:
            messagebox.showerror("Pipeline failed", str(exc))
            log(f"ERROR: {exc}")

    tk.Radiobutton(mode_frame, text="Genome FASTA + GFF3", variable=mode_var,
                   value="genome_gff", command=toggle_mode).grid(row=0, column=0, sticky="w")
    tk.Radiobutton(mode_frame, text="CDS FASTA only", variable=mode_var,
                   value="cds_only", command=toggle_mode).grid(row=0, column=1, sticky="w")

    tk.Button(root, text="🚀 Design Primers", command=run_pipeline,
              bg="lightgreen", font=("Arial", 14, "bold"), height=2).grid(
        row=4, column=0, columnspan=6, pady=12)

    toggle_mode()
    for i in range(6):
        root.columnconfigure(i, weight=1)
    root.mainloop()


if __name__ == "__main__":
    main()
