# Improvement Plan & Roadmap

This document tracks what was fixed/added and what is recommended next, ordered
by impact on the project's two priorities: **accuracy** and **ease of use**.

---

## ✅ Done in this iteration (latest)

### Accuracy — off-target mismatch reporting (NEW)
- **Each predicted amplicon now carries its mismatch count.** `in_silico_pcr()`
  amplicons are `(chrom, start, end, size, mismatches)`; the intended product is
  0 mismatches and off-targets are ≥1 in mismatch-tolerant mode.
  `specificity_label()` surfaces the **nearest off-target's mismatch count**, so
  the CSV `Specificity Check` now distinguishes "non-specific, but the nearest
  off-target carries 3 mismatches" from a perfect-match off-target — the former
  is far more likely to be fine at the bench. Regression-guarded
  (`test_specificity_label_reports_offtarget_mismatches`,
  `test_in_silico_pcr_amplicon_carries_mismatch_count`).

### Accuracy bug fix — amplicon enumeration early-break (NEW)
- The amplicon search sorted reverse binding sites by position and `break`-ed on
  the first product over `max_product`. When the two primers differ in length a
  later, longer-footprint site could still be in-window, so a valid (often
  larger) off-target amplicon could be **silently skipped**. The break now uses
  the *minimum* reverse footprint length as a correct lower bound, so no
  in-window amplicon is missed while the prune is preserved. Regression-guarded
  (`test_in_silico_pcr_unequal_primer_lengths_large_product`); benchmarked — no
  speed regression (specificity path unchanged within noise).

### Ease of use + accuracy — ranked alternate primer pairs (NEW)
- **Report the top *N* primer pairs per target, not just the single best.**
  `PrimerParams.num_return` is now honoured end-to-end: `design_primer_candidates()`
  returns ranked results (rank 0 = best), `design_for_gene()` takes
  `num_candidates`, and a `Rank` column was added to the CSV. Exposed as
  `--num-return` (CLI) and *Primers/target* (GUI). If the best pair fails at the
  bench the user has scored fallbacks without re-running the pipeline.
  `design_primers_for_sequence()` is kept as a back-compatible single-result
  wrapper. Regression-guarded (`test_design_primer_candidates_ranked_alternates`,
  `test_design_primers_for_sequence_is_first_candidate`,
  `test_design_for_gene_num_candidates_multiplies_rows`) and benchmarked (new
  "Ranked alternates" section in `benchmark.py`: 3 pairs ≈ 1.5× the time of 1 for
  3× the output, since it is a single Primer3 call).

---

## ✅ Done in previous iteration

### Accuracy — 3′-anchored, mismatch-tolerant specificity
- **Off-targets bind with mismatches; exact matching under-reports them.** The
  in-silico PCR check now offers a 3′-anchored, mismatch-tolerant search:
  `in_silico_pcr(..., seed_len=12, max_mismatches=2)`. A binding site counts
  when its **3′-most `seed_len` bases match exactly** (the seed) and the whole
  footprint differs from the primer by at most `max_mismatches` bases. This
  mirrors the biology — a primer extends only when its 3′ end is matched, while
  5′ mismatches are tolerated — so it surfaces off-targets the exact search
  silently missed. `seed_len=0` (default) keeps the exact behaviour bit-for-bit.
  - Implementation seeds on the exact 3′ k-mer with `str.find` and only scores
    the few candidate footprints, so it stays within **~9 % of the exact path**
    (14.4 vs 13.2 ms/pair at *E. coli* scale) and still ~2× faster than the
    original reload-per-pair code. See `benchmark.py`.
  - Exposed in the CLI (`--seed-len` / `--max-mismatches`) and GUI (*3′ seed* /
    *Max mismatches* fields), and recorded in the CSV parameter banner.
  - Regression-guarded: `test_in_silico_pcr_mismatch_offtarget_detected`
    (5′-mismatched off-target found), `test_in_silico_pcr_3prime_mismatch_not_extended`
    (3′-seed mismatch correctly rejected), `test_in_silico_pcr_seed_len_zero_is_exact`
    (exact-mode equivalence).

### Primer placement control
- **Choose where each primer lands relative to the gene.** Via Primer3's
  `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`, the forward and reverse primers can be
  constrained independently to the **upstream flank**, **inside the gene**, or
  the **downstream flank**. Modes: `internal` (default, back-compatible),
  `flanking` (upstream→downstream, e.g. knockout verification), `custom`
  (any region pair), and `all` (every valid permutation, one output row each).
  Implemented in `build_gene_template()`, `placement_combos()` and
  `design_for_gene()`; exposed in both the GUI (dropdown) and CLI
  (`--placement` / `--fwd-region` / `--rev-region`). Regression-tested
  (`test_design_flanking_places_primers_in_flanks`, `test_design_all_permutations`)
  and benchmarked (placement section in `benchmark.py`).
- Template coordinates are preserved for region constraints (`_p3_template`
  replaces non-ACGT with `N` instead of deleting it, so offsets stay valid).

### Accuracy
- **Two-orientation in-silico PCR.** The original specificity check only looked
  for the forward primer on the top strand plus the reverse primer's
  reverse-complement on the top strand. It silently missed off-target amplicons
  where the primers' roles are swapped (reverse primer acting as the left
  primer) and self-amplicons. `in_silico_pcr()` now enumerates every binding
  site of **both** primers in **both** senses and counts every product that
  falls inside the size window. Regression-guarded by
  `tests/test_in_silico_pcr_reverse_orientation`.
- **Internally consistent Tm.** Reported Tm now uses the *same* salt / dNTP /
  DNA concentrations that Primer3 used during design (`ThermoParams`), instead
  of Primer3's library defaults. Previously the reported Tm could differ from
  the Tm the design was optimised against.
- **Secondary-structure metrics.** Each result now reports hairpin Tm,
  self-dimer (homodimer) Tm for both primers, and primer-pair hetero-dimer Tm
  (`primer3.calc_hairpin/homodimer/heterodimer`). These flag primers that fold
  or dimerise instead of binding template — a common cause of failed PCR.
- **GC-clamp** constraint exposed (default 1) to favour a 3′ G/C and improve
  priming efficiency.
- **Honest failure reporting.** Templates shorter than the primer or the minimum
  product size now return an explanatory status instead of a generic
  "No suitable primers", and Primer3's `PRIMER_PAIR_EXPLAIN` is surfaced.

### Correctness bugs fixed
- **Gene-filter placeholder bug.** The default example text
  (`# sulA, opgH, galU`) was parsed as real gene names because `#` comment lines
  were not ignored — so leaving the box "empty" silently filtered to those
  genes. `parse_target_names()` now ignores `#` lines and blanks, de-duplicates,
  and lower-cases. Regression-guarded.
- **Specificity genome reload.** The genome FASTA was re-parsed **and re-cleaned
  from disk for every primer pair**. It is now loaded once and cleaned once
  (`prepare_genome()`), giving a ~2x end-to-end speed-up at *E. coli* scale
  while doing strictly more work (two orientations). See `benchmark.py`.
- **Bare `except:` everywhere** replaced with narrow handling + logging, so real
  errors are no longer swallowed.

### Ease of use
- **Importable / testable core.** All logic moved to `primer_design.py` with no
  Tkinter dependency. The GUI (`enhanced_primer_gui.py`) no longer launches on
  import — it only runs under `__main__` — so the logic can be scripted and
  tested headlessly.
- **Headless CLI** (`primer_cli.py`) for batch runs and automation, sharing
  100% of its logic with the GUI.
- **Parameter validation** with clear messages before a run starts.
- **Reproducible env**: `requirements.txt` + `environment.yml` (the CI workflow
  referenced an `environment.yml` that did not exist — CI was broken).
- **Test suite** (`tests/`, 17 tests) and a **benchmark harness**
  (`benchmark.py`) so every future edit can be validated for correctness and
  speed.

---

## 🔜 Recommended next (high value)

1. **Gold-standard specificity via BLAST.** The 3′-anchored mismatch search
   (done) plus per-amplicon mismatch reporting (done — see above) cover the
   common case. For publication-grade off-target prediction, integrate NCBI
   Primer-BLAST or a local BLAST/`isPcr` so indels and degenerate sites are
   handled too.
2. **Faster specificity at genome scale.** The current search is `str.find`
   over each contig per pair (O(genome × pairs)). Build a k-mer index or use a
   suffix automaton / Aho-Corasick over all primers at once for large genomes
   or whole-genome primer panels.
3. **Run design off the main thread.** For large gene sets the Tkinter GUI
   freezes during a run; move the pipeline to a worker thread with a progress
   queue.
4. **Rank/score the alternates.** Ranked alternates are now emitted (done — see
   above); a natural follow-up is a composite quality score per pair (penalising
   hairpin/dimer Tm near the annealing temp, GC/Tm imbalance) to re-rank or flag
   risky candidates beyond Primer3's default ordering.

## 💡 Nice to have
- **Expose warning thresholds** (`WARN_MAX_TM_DIFF`, `WARN_HAIRPIN_TM`,
  `WARN_DIMER_TM`) as CLI flags / GUI fields so users can tune the sensitivity
  of the `Warnings` column to their assay, and optionally let a warning demote a
  candidate so a cleaner alternate is chosen (ties into auto-pick above).
- Save/load parameter presets (JSON) from the GUI.
- Allow primers to be placed in flanking regions (knockout-verification
  primers), not just inside the gene; the flank size is already extracted.
- Export GenBank/BED of amplicon locations.
- Amplicon Tm / size-distribution plots.
- Package as `pip install`-able with a console entry point.

---

## How to validate changes

```bash
pip install -r requirements.txt
pytest -q                 # correctness
python benchmark.py       # speed (writes benchmark_results.md)
flake8 . --select=E9,F63,F7,F82   # CI-blocking lint
```
