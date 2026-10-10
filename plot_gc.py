
"""Plot sliding-window GC content from a SeqToolkit CSV report."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def plot_gc_profile(csv_path, output_path):
    """Plot GC percentage against the midpoint of each window."""
    csv_path = Path(csv_path)

    positions = []
    gc_values = []
    current_header = None

    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)

        required = {"header", "start", "end", "GC_percent"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("CSV is missing required columns.")

        for row in reader:
            header = row["header"]

            # This first version plots one sequence at a time.
            if current_header is None:
                current_header = header

            if header != current_header:
                continue

            positions.append(
                (int(row["start"]) + int(row["end"])) / 2
            )
            gc_values.append(float(row["GC_percent"]))

    if not positions:
        raise ValueError("No GC-window data found in the CSV.")

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(positions, gc_values, marker="o", linewidth=2)
    ax.set_title(f"GC Content Profile: {current_header}")
    ax.set_xlabel("Sequence position (bp; window midpoint)")
    ax.set_ylabel("GC content (%)")
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"GC profile saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot a SeqToolkit sliding-window GC report."
    )
    parser.add_argument(
        "--input",
        default="output/gc_windows.csv",
        help="Path to the GC-window CSV",
    )
    parser.add_argument(
        "--output",
        default="output/gc_profile.png",
        help="Path for the plot image",
    )
    args = parser.parse_args()

    try:
        plot_gc_profile(args.input, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
