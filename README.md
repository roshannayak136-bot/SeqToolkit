
# SeqToolkit 🧬

A Python command-line toolkit for basic DNA sequence analysis
and sliding-window GC-content visualization.

## Features

- Parse multi-record FASTA files
- Calculate sequence length and nucleotide counts
- Calculate GC percentage
- Count ambiguous DNA symbols
- Generate reverse complements
- Calculate sliding-window GC content
- Export analysis results to CSV
- Plot GC-content profiles
- Run automated unit tests

## Requirements

- Python 3.10+
- matplotlib for visualization

## Installation

```bash
git clone https://github.com/roshannayak136-bot/SeqToolkit.git
cd SeqToolkit
python -m venv .venv
```

On Windows PowerShell:

```powershell
.\.venv\Scripts\python.exe -m pip install matplotlib
```

## Run the analysis

```bash
python seqtoolkit.py data/sample.fasta --window-size 10 --step-size 5
```

Generate a GC profile:

```bash
python plot_gc.py
```

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Scientific limitations

The included FASTA sequences are synthetic demonstration data.
GC calculations include ambiguous bases in the denominator,
so the reported percentage is not an estimate of true GC content
when ambiguous bases are present.

This educational toolkit is not validated for clinical use and
does not perform sequence alignment, variant calling, or
clinical interpretation.

## Author

MSc Biotechnology portfolio project by Roshan Nayak.

