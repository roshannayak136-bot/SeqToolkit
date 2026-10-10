
"""Basic quality-control checks for DNA FASTA records."""

import argparse
import csv
from pathlib import Path

from seqtoolkit import parse_fasta, VALID_DNA

UNAMBIGUOUS_DNA = set("ACGT")


def quality_check(header, sequence, min_length=10, max_length=100000):
    """Return QC metrics and warnings for one DNA sequence."""
    sequence = sequence.upper()

    if not sequence:
        raise ValueError(f"{header}: sequence is empty")

    invalid = sorted(set(sequence) - VALID_DNA)
    if invalid:
        raise ValueError(
            f"{header}: invalid DNA symbols: {', '.join(invalid)}"
        )

    length = len(sequence)
    ambiguous = sum(base not in UNAMBIGUOUS_DNA for base in sequence)
    gc_bases = sequence.count("G") + sequence.count("C")
    canonical_bases = length - ambiguous

    # GC among A/C/G/T bases only; undefined if none are present.
    gc_percent = (
        round(gc_bases / canonical_bases * 100, 2)
        if canonical_bases else None
    )

    warnings = []

    if length < min_length:
        warnings.append("Sequence shorter than minimum length")

    if length > max_length:
        warnings.append("Sequence longer than maximum length")

    if ambiguous:
        warnings.append(f"{ambiguous} ambiguous base(s)")

    status = "WARNING" if warnings else "PASS"

    return {
        "header": header,
        "length": length,
        "ambiguous_count": ambiguous,
        "ambiguous_percent": round(ambiguous / length * 100, 2),
        "GC_percent_unambiguous": gc_percent,
        "status": status,
        "warnings": "; ".join(warnings) if warnings else "None",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run basic quality checks on DNA FASTA records."
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument(
        "--output", default="output/qc_report.csv",
        help="Output CSV path"
    )
    parser.add_argument(
        "--min-length", type=int, default=10,
        help="Minimum expected sequence length"
    )
    parser.add_argument(
        "--max-length", type=int, default=100000,
        help="Maximum expected sequence length"
    )
    args = parser.parse_args()

    if args.min_length < 1:
        parser.error("--min-length must be at least 1")
    if args.max_length < args.min_length:
        parser.error("--max-length must be >= --min-length")

    try:
        records = parse_fasta(args.input)
        results = [
            quality_check(
                header, sequence,
                min_length=args.min_length,
                max_length=args.max_length,
            )
            for header, sequence in records
        ]

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

    except (OSError, ValueError) as error:
        parser.error(str(error))

    for result in results:
        print(
            f'{result["status"]:7} | {result["header"]} | '
            f'length={result["length"]} | '
            f'ambiguous={result["ambiguous_count"]} | '
            f'GC={result["GC_percent_unambiguous"]}% | '
            f'{result["warnings"]}'
        )

    print(f"\nQC report saved to: {output_path}")


if __name__ == "__main__":
    main()
