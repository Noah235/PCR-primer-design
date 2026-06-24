# Improvement Plan & Roadmap

This document tracks what was fixed/added and what is recommended next, ordered
by impact on the project's two priorities: **accuracy** and **ease of use**.

---

## ✅ Done in this iteration

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

1. **Mismatch-tolerant / 3′-anchored specificity.** Real off-targets bind with
   mismatches. Exact matching under-reports non-specificity. Options, in order
   of effort:
   - Seed on the 3′-most ~12 nt (most important for extension), then score
     mismatches in the remainder.
   - Integrate NCBI Primer-BLAST or a local BLAST/`isPcr` for gold-standard
     specificity.
2. **Faster specificity at genome scale.** The current search is `str.find`
   over each contig per pair (O(genome × pairs)). Build a k-mer index or use a
   suffix automaton / Aho-Corasick over all primers at once for large genomes
   or whole-genome primer panels.
3. **Run design off the main thread.** For large gene sets the Tkinter GUI
   freezes during a run; move the pipeline to a worker thread with a progress
   queue.
4. **Multiple primer candidates per gene.** `PrimerParams.num_return` is already
   plumbed through; expose ranked alternates (PRIMER_*_1, _2 …) in the output.

## 💡 Nice to have
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
