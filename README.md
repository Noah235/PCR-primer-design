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
- **3′-anchored, mismatch-tolerant specificity** (optional) — models the biology
  that a primer only extends when its 3′ end matches: the 3′ seed must match
  exactly while a budget of 5′ mismatches is allowed, surfacing off-targets that
  exact matching silently misses. Enabled with `--seed-len`/`--max-mismatches`
  (CLI) or the *3′ seed* / *Max mismatches* fields (GUI).
- **Secondary-structure reporting** per primer: hairpin Tm, self-dimer Tm, and
  primer-pair hetero-dimer Tm — so you can spot primers likely to fail.
- **Internally consistent Tm** (reported under the same salt conditions used for
  design).
- **Primer placement control** — choose where each primer lands relative to the
  gene: both inside the gene (default), forward in the upstream flank + reverse
  in the downstream flank (*flanking*, e.g. knockout verification), any custom
  upstream/internal/downstream combination, or **all permutations** at once
  (one row per placement).
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

# Primer placement: amplify each gene from upstream flank to downstream flank
python primer_cli.py genome --genome genome.fasta --gff annotation.gff3 \
    --placement flanking --flank 250 -o flanking_primers.csv

# Every placement permutation (6 rows per gene), or a custom combination
python primer_cli.py genome ... --placement all
python primer_cli.py genome ... --placement custom --fwd-region upstream --rev-region internal

# Mismatch-tolerant specificity: 12 nt exact 3' seed, up to 2 mismatches in the 5' tail
python primer_cli.py genome --genome genome.fasta --gff annotation.gff3 \
    --specificity --seed-len 12 --max-mismatches 2 -o primers.csv

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

# Mismatch-tolerant: catch off-targets with 5' mismatches (3' seed must match)
spec = pd.in_silico_pcr(result["forward"], result["reverse"], genome,
                        seed_len=12, max_mismatches=2)
```

---

## Output (CSV columns)

Gene name · Placement · Forward/Reverse primer · Tm · GC% · Product length ·
Hairpin Tm (F/R) · Self-dimer Tm (F/R) · Hetero-dimer Tm ·
Specificity check · Status.

The **Placement** column records where the pair was designed (e.g.
`internal->internal`, `upstream->downstream`). With `--placement all` you get one
row per permutation per gene.

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

- The in-silico specificity check defaults to **exact** matching (fast screening).
  For higher accuracy enable the **3′-anchored mismatch-tolerant** mode
  (`--seed-len 12 --max-mismatches 2`), which finds off-targets carrying 5′
  mismatches. Experimental validation is still recommended; for gold-standard
  results integrate Primer-BLAST / local BLAST (see `IMPROVEMENTS.md`).
- Specificity scales as O(genome × primer pairs); for whole-genome panels a
  k-mer index is the recommended next step.

## How to extend

See **`IMPROVEMENTS.md`** for a prioritised roadmap.
