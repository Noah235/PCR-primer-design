# Improvement Plan & Roadmap

This document tracks what was fixed/added and what is recommended next, ordered
by impact on the project's two priorities: **accuracy** and **ease of use**.

---

## ✅ Done in this iteration (latest)

### Ease of use + accuracy — `check` mode for primers you already have (NEW)
- **The tool could only score primers it *designed*.** A bench scientist with an
  existing primer pair (from a paper, a collaborator, or an earlier run) had no
  way to run the same QC or the accurate in-silico-PCR specificity check without
  reverse-engineering a template to feed the design pipeline — a real, common
  workflow was simply impossible.
- **Fix / feature:** a new `analyze_primer_pair()` computes the full diagnostic
  set for an arbitrary pair — each primer's Tm and GC%, hairpin / self-dimer /
  hetero-dimer Tm, the composite quality score and the at-a-glance warnings —
  reusing the exact same functions the design pipeline uses (so a checked primer
  reports bit-identical Tm/GC to a designed one). A new `read_primer_pairs()`
  loads a batch file of `name forward reverse` lines (comma/tab/space separated,
  `#` comments and blanks ignored, 2-field lines auto-named). Both are wired into
  a new **`check` CLI subcommand**: pass `--forward/--reverse` or `--primers
  FILE`, optionally a `--genome` to run the two-orientation (optionally
  mismatch-tolerant, `--seed-len/--max-mismatches`) specificity search, and
  optionally `--bed` to export every predicted amplicon. With a genome the
  predicted **product size** is filled from the in-silico PCR result and a single
  predicted amplicon is reported as *Specific*; without a genome it is pure
  primer QC (still fully useful). Output uses the shared CSV column layout.
- **Regression-guarded** (6 new tests: metric reporting + emittable row,
  designed-vs-checked Tm/GC identity, empty/too-short input flagged, batch-file
  parsing across comma/tab/space + 2-field + malformed-skip, `check` CLI
  end-to-end reporting a specific pair as one amplicon, and `check` QC-only with
  no genome). **Benchmarked** (new "Primer check" section): QC alone is
  ~0.5 ms/pair (a few secondary-structure calls, no Primer3 design) and adding
  the genome specificity check costs the same ~14 ms/pair two-orientation search
  the design pipeline already uses — the new front-end adds no algorithmic cost.
  Verified end-to-end on a synthetic genome carrying a **duplicated** locus: the
  pair is correctly flagged *Non-specific (3 amplicons)* and every predicted
  amplicon (both real copies plus the spurious cross-copy product) is written to
  BED, while a deliberately too-short pair fails gracefully with an explanatory
  status instead of crashing the batch.

### Accuracy + ease of use — amplicon location BED export (NEW)
- **You could see *that* a pair was non-specific, but not *where*.** The
  specificity search already computes every predicted amplicon's genome
  coordinates, but they were collapsed into a one-line label (`Non-specific
  (3 amplicons)`) and then discarded — so a user had no way to inspect the
  off-targets, judge whether they matter, or design around them.
- **Fix / feature:** a new `amplicon_bed_rows()` turns the amplicon tuples into
  BED6 rows and `write_bed()` writes a browser-ready, coordinate-sorted BED with
  a UCSC track header. Each feature is named
  `{gene}_rank{N}_{placement}_{ontarget|offtarget}_mm{k}` and shaded by mismatch
  count (score `1000 − 250·mm`), so a genome browser (IGV/JBrowse/UCSC) shows the
  intended product **and** every off-target at a glance, darker = more likely to
  prime. On/off-target is decided by overlap with the gene's own locus, so the
  intended amplicon is labelled automatically. Exposed as `--bed FILE` (CLI,
  guarded to require `--specificity` and failing fast before the genome loads)
  and an *Amplicon BED (optional)* field in the GUI.
- **Regression-guarded** (7 new tests: on/off-target classification incl.
  wrong-contig, no-target → all off-target, 1-based rank + placement in the
  feature name, empty-amplicon → no rows, `write_bed` sorts + writes the header
  and BED6 layout, empty export still writes a valid header-only file, plus a
  `specificity_label` singular-amplicon guard). **Benchmarked** (new "Amplicon
  BED export" section: ~1.9 µs/amplicon — a pure in-memory transform of tuples
  the specificity search already produced, so `--bed` adds no meaningful cost).
  Verified end-to-end on a synthetic genome carrying a duplicated gene: the
  duplicate is correctly emitted as an `offtarget` feature at its coordinates
  while the real locus is `ontarget`.
- **Also fixed** a dead computed variable (`note`) in `specificity_label()` that
  was assigned but never used.

### Bug fix (ease of use + accuracy) — tolerant FASTA loading (NEW)
- **Real genome / CDS FASTAs crashed the loader on current Biopython.**
  Biopython **>= 1.85** made the default `"fasta"` parser *strict*: it now
  raises on any file that has blank or `;`-comment lines **before the first
  record** (and on `;` comment lines anywhere). Such lines are common in
  real-world FASTAs (and the Pearson/BLAST FASTA dialects allow them), so
  `load_genome()` / `load_cds_sequences()` — which both called
  `SeqIO.parse(path, "fasta")` — failed with an opaque multi-paragraph
  Biopython traceback on input that loaded fine on older Biopython. The whole
  genome+GFF3 pipeline (CLI **and** GUI) was unusable on those files.
- **Fix:** a new `read_fasta_records()` tries the fast strict parser first and,
  only when it raises, falls back to stripping the offending leading/comment
  lines (`_strip_leading_fasta_comments()`) and re-parsing the cleaned text.
  This is version-agnostic (no dependency on the new `fasta-pearson` /
  `fasta-blast` format names, which don't exist on Biopython 1.80–1.84) and
  both loaders now route through it. Empty/headerless files now surface the
  clear `No sequences found in FASTA file` message instead of the strict
  parser's traceback.
- **Regression-guarded** (5 new tests: leading-comment genome loads, clean
  multi-record file still loads via the fast path, empty file → clear error,
  CDS leading-blank + inline `;` comment loads, `_strip_leading_fasta_comments`
  unit). **Benchmarked** (new "FASTA loading" section): the strict fast path is
  unchanged for clean files (~1.5 ms / 200 records) and the fallback (~3.2 ms)
  only runs for files that would otherwise crash — **no overhead on the common
  path**. Verified end-to-end: the full `genome` CLI pipeline (specificity +
  ranked alternates) now completes on a FASTA with a leading comment block that
  previously aborted before any primer was designed.

### Accuracy + ease of use — composite quality score & quality re-ranking
- **One 0–100 score per pair, and the option to rank by it.** `quality_score()`
  folds the numbers a bench scientist would otherwise weigh by eye — each
  primer's Tm deviation from `opt_tm`, the pair Tm difference, GC distance from
  50 %, and how stable the hairpin / self-dimer / hetero-dimer structures are
  *relative to the annealing temperature* — into a single figure where higher is
  better (structures that melt out well before annealing are not penalised).
  Surfaced as a new `Quality Score` CSV column and computed for **every**
  candidate during design (no extra Primer3 calls).
- **`rank_by="quality"`** (CLI `--rank-by-quality`, GUI *Rank by quality score*)
  re-orders the alternates by this score so a structurally cleaner pair can be
  promoted to rank 1 ahead of Primer3's default ordering — which does **not**
  factor in the hairpin/dimer Tm values we compute. Default ordering is
  unchanged (`rank_by="primer3"`), so existing output is bit-for-bit stable.
  Threaded through `design_primer_candidates`, `design_primers_for_sequence`,
  `design_for_gene`, the CLI and the GUI. Regression-guarded (8 new tests:
  clean-pair=100, monotonic penalties, `None` for failure rows, clamp-at-0,
  end-to-end score in the CSV row, quality ordering is a sorted permutation,
  rank-0 carries the best score, invalid `rank_by` raises) and benchmarked (new
  "Quality re-ranking" section: re-ranking is an in-memory sort, **0 % runtime
  overhead** vs plain ranked alternates — 10.9 vs 11.0 ms/template).

### Ease of use — GUI dead-widget fix (NEW)
- The GUI created the *Candidates* / *Primers per target* number entry **twice**
  (`num_return_e` was assigned over itself); the first widget was orphaned and
  never read. Removed the duplicate so the visible field is the one actually used.

### Accuracy — off-target mismatch reporting
- **Each predicted amplicon now carries its mismatch count.** `in_silico_pcr()`
  amplicons are `(chrom, start, end, size, mismatches)`; the intended product is
  0 mismatches and off-targets are ≥1 in mismatch-tolerant mode.
  `specificity_label()` surfaces the **nearest off-target's mismatch count**, so
  the CSV `Specificity Check` now distinguishes "non-specific, but the nearest
  off-target carries 3 mismatches" from a perfect-match off-target — the former
  is far more likely to be fine at the bench. Regression-guarded
  (`test_specificity_label_reports_offtarget_mismatches`,
  `test_in_silico_pcr_amplicon_carries_mismatch_count`).

### Accuracy bug fix — amplicon enumeration early-break
- The amplicon search sorted reverse binding sites by position and `break`-ed on
  the first product over `max_product`. When the two primers differ in length a
  later, longer-footprint site could still be in-window, so a valid (often
  larger) off-target amplicon could be **silently skipped**. The break now uses
  the *minimum* reverse footprint length as a correct lower bound, so no
  in-window amplicon is missed while the prune is preserved. Regression-guarded
  (`test_in_silico_pcr_unequal_primer_lengths_large_product`); benchmarked — no
  speed regression (specificity path unchanged within noise).

### Ease of use + accuracy — ranked alternate primer pairs
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
4. **Rank/score the alternates.** ✅ *Done* — a composite 0–100 `quality_score()`
   (Tm match/balance, GC, hairpin/dimer Tm relative to the annealing temp) is now
   computed per pair and reported in a `Quality Score` column, and
   `rank_by="quality"` (`--rank-by-quality`) re-ranks the alternates by it beyond
   Primer3's default ordering. Natural follow-ups: expose the penalty weights /
   structure margin as tunables, and let a low score demote or flag a candidate
   automatically (ties into the warning thresholds below).

## 💡 Nice to have
- **Expose warning thresholds** (`WARN_MAX_TM_DIFF`, `WARN_HAIRPIN_TM`,
  `WARN_DIMER_TM`) as CLI flags / GUI fields so users can tune the sensitivity
  of the `Warnings` column to their assay, and optionally let a warning demote a
  candidate so a cleaner alternate is chosen (ties into auto-pick above).
- Save/load parameter presets (JSON) from the GUI.
- Allow primers to be placed in flanking regions (knockout-verification
  primers), not just inside the gene; the flank size is already extracted.
- Export GenBank/BED of amplicon locations. ✅ *BED done* (`--bed` / GUI
  *Amplicon BED*; on/off-target-classified, mismatch-shaded, browser-ready).
  Follow-ups: GenBank export, and BED for CDS mode (needs a coordinate frame).
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
