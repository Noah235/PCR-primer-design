# PCR Primer Design with Specificity Testing

A Python tool for designing PCR primers from a reference genome + annotations
(or a CDS FASTA), with in-silico specificity testing. Usable as a **Tkinter
GUI** or a **headless CLI**, with all logic in an importable, tested core module.

---

## Features

- Two input modes:
  - **Genome FASTA + GFF3** — extract each gene (with optional flanks) and design primers.
  - **CDS FASTA only** — design directly from coding sequences.
- Primer design via Primer3 with tunable size, Tm, GC%, product size and GC-clamp.
- **Accurate, two-orientation in-silico PCR** specificity check that counts every
  predicted amplicon across the genome (catches off-targets the naive
  single-orientation search misses).
- **Secondary-structure reporting** per primer: hairpin Tm, self-dimer Tm, and
  primer-pair hetero-dimer Tm — so you can spot primers likely to fail.
- **Internally consistent Tm** (reported under the same salt conditions used for
  design).
- Case-insensitive gene/locus-tag filtering (`#` comment lines are ignored).
- CSV output + extracted-sequence FASTA.

---

## Installation

```bash
pip install -r requirements.txt
# or, with conda:
conda env create -f environment.yml && conda activate pcr-primer-design
```

Requires Python 3.9+, `biopython` and `primer3-py`. The GUI additionally needs
`tkinter` (bundled with most Python installs; on Debian/Ubuntu:
`sudo apt install python3-tk`). The CLI and core library do **not** need tkinter.

---

## Usage

### GUI

```bash
python enhanced_primer_gui.py
```

Pick an input mode, select files, set parameters, optionally enable specificity
testing, and click **Design Primers**.

### CLI (headless / batch)

```bash
# Genome + GFF3, three genes, with specificity testing
python primer_cli.py genome \
    --genome genome.fasta --gff annotation.gff3 \
    --genes sulA,opgH,galU --specificity -o primers.csv

# CDS FASTA, all records
python primer_cli.py cds --cds cds.fasta -o primers.csv

python primer_cli.py genome --help    # full parameter list
```

### As a library

```python
import primer_design as pd

params = pd.PrimerParams(opt_tm=60, product_min=100, product_max=600)
result = pd.design_primers_for_sequence("myGene", template_seq, params)

genome = pd.prepare_genome(pd.load_genome("genome.fasta"))
spec = pd.in_silico_pcr(result["forward"], result["reverse"], genome)
print(pd.specificity_label(spec))
```

---

## Output (CSV columns)

Gene name · Forward/Reverse primer · Tm · GC% · Product length ·
Hairpin Tm (F/R) · Self-dimer Tm (F/R) · Hetero-dimer Tm ·
Specificity check · Status.

In genome mode an additional `*_extracted_sequences.fasta` is written with each
gene's 5′ flank, 3′ flank and coding region.

---

## Project layout

| File | Purpose |
| --- | --- |
| `primer_design.py` | Core logic (no GUI dependency) — import & test this |
| `primer_cli.py` | Command-line front-end |
| `enhanced_primer_gui.py` | Tkinter GUI front-end |
| `tests/` | pytest correctness suite |
| `benchmark.py` | Speed benchmark (writes `benchmark_results.md`) |
| `IMPROVEMENTS.md` | Roadmap and changelog |

---

## Development

```bash
pytest -q                          # run tests
python benchmark.py                # measure performance
flake8 . --select=E9,F63,F7,F82    # CI-blocking lint
```

---

## Notes & limitations

- The in-silico specificity check uses **exact** matching (no mismatches). It is
  fast and good for screening, but experimental validation is still recommended.
  See `IMPROVEMENTS.md` for the plan to add mismatch-tolerant / BLAST-based
  specificity.
- Specificity scales as O(genome × primer pairs); for whole-genome panels a
  k-mer index is the recommended next step.

## How to extend

See **`IMPROVEMENTS.md`** for a prioritised roadmap.
