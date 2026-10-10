
"""SeqToolkit: a small, reproducible DNA sequence analysis tool."""

import argparse
import csv
from pathlib import Path


# IUPAC DNA alphabet: A, C, G, T plus ambiguity codes.
VALID_DNA = set("ACGTRYSWKMBDHVN")

# Each ambiguity symbol is paired with its complement.
COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


def parse_fasta(path):
    """Read FASTA records and return a list of (header, sequence) tuples."""
    records = []
    header = None
    sequence_parts = []

    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_parts).upper()
                    if not sequence:
                        raise ValueError(
                            f"Empty sequence for FASTA record: {header}"
                        )
                    records.append((header, sequence))

                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Missing FASTA header at line {line_number}"
                    )
                sequence_parts = []

            else:
                if header is None:
                    raise ValueError(
                        f"Sequence found before FASTA header "
                        f"at line {line_number}"
                    )

                # Permit whitespace within sequence lines.
                sequence_parts.append("".join(line.split()))

    if header is not None:
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record: {header}")
        records.append((header, sequence))

    if not records:
        raise ValueError("No FASTA records found.")

    return records


def analyze_sequence(header, sequence):
    """Calculate composition and reverse complement for one DNA sequence."""
    sequence = sequence.upper()

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    invalid = sorted(set(sequence) - VALID_DNA)
    if invalid:
        raise ValueError(
            f"Invalid DNA symbols in {header}: {', '.join(invalid)}"
        )

    length = len(sequence)
    counts = {base: sequence.count(base) for base in "ACGT"}

    # GC% uses the full sequence length as denominator.
    # Ambiguous symbols remain in the denominator.
    gc_percent = (counts["G"] + counts["C"]) / length * 100

    return {
        "header": header,
        "length": length,
        "A_count": counts["A"],
        "C_count": counts["C"],
        "G_count": counts["G"],
        "T_count": counts["T"],
        "ambiguous_count": sum(
            base not in "ACGT" for base in sequence
        ),
        "GC_percent": round(gc_percent, 2),
        "reverse_complement": sequence.translate(COMPLEMENT)[::-1],
    }


def write_csv(results, output_path):
    """Write analysis results to a CSV file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(results[0].keys())

    with output_path.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

def calculate_gc_windows(sequence, window_size=10, step_size=5):
    """Calculate GC percentage across overlapping DNA windows."""
    sequence = sequence.upper()

    if window_size <= 0 or step_size <= 0:
        raise ValueError("Window size and step size must be positive.")

    invalid = sorted(set(sequence) - VALID_DNA)
    if invalid:
        raise ValueError(
            f"Invalid DNA symbols: {', '.join(invalid)}"
        )

    results = []

    for start in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[start:start + window_size]
        gc_count = window.count("G") + window.count("C")

        results.append({
            "start": start + 1,
            "end": start + window_size,
            "GC_percent": round(gc_count / window_size * 100, 2),
        })

    return results



def main():
    parser = argparse.ArgumentParser(
        description="Analyze DNA sequences from a FASTA file."
    )
    parser.add_argument(
        "input",
        help="Path to the input FASTA file",
    )
    parser.add_argument(
        "--output",
        default="output/summary.csv",
        help="Path for the sequence summary CSV",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=10,
        help="Sliding-window length in bases (default: 10)",
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=5,
        help="Distance between window starts (default: 5)",
    )
    args = parser.parse_args()

    try:
        records = parse_fasta(args.input)

        results = [
            analyze_sequence(header, sequence)
            for header, sequence in records
        ]
        write_csv(results, args.output)

        window_results = []

        for header, sequence in records:
            windows = calculate_gc_windows(
                sequence,
                window_size=args.window_size,
                step_size=args.step_size,
            )

            for window in windows:
                window_results.append({
                    "header": header,
                    **window,
                })

        window_output = Path(args.output).with_name(
            "gc_windows.csv"
        )
        if window_results:
            write_csv(window_results, window_output)

    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"Analyzed {len(results)} DNA sequence(s).")
    print(f"Sequence summary saved to: {args.output}")

    if window_results:
        print(f"GC window report saved to: {window_output}")
    else:
        print("No complete GC windows found.")


if __name__ == "__main__":
    main()



