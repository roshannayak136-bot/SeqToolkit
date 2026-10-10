
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def plot_gc_profiles(csv_path, output_dir):
    """Create a separate GC-content plot for each sequence in a CSV."""
    csv_path = Path(csv_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    grouped = defaultdict(lambda: {"positions": [], "gc_values": []})

    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)

        required = {"header", "start", "end", "GC_percent"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("CSV is missing required columns.")

        for row in reader:
            header = row["header"]

            try:
                start = int(row["start"])
                end = int(row["end"])
                gc = float(row["GC_percent"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid GC-window data for sequence {header!r}."
                ) from error

            if start < 0 or end <= start or not 0 <= gc <= 100:
                raise ValueError(
                    f"Out-of-range GC-window data for sequence {header!r}."
                )

            grouped[header]["positions"].append((start + end) / 2)
            grouped[header]["gc_values"].append(gc)

    if not grouped:
        raise ValueError("No GC-window data found in the CSV.")

    output_paths = []

    for index, (header, values) in enumerate(grouped.items(), start=1):
        # Make a filesystem-friendly filename.
        slug = re.sub(r"[^A-Za-z0-9_-]+", "_", header).strip("_")
        slug = slug[:60] or "sequence"
        output_path = output_dir / f"gc_profile_{index:03d}_{slug}.png"

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(
            values["positions"],
            values["gc_values"],
            marker="o",
            linewidth=1.8,
        )
        ax.set_title(f"GC Content Profile: {header}")
        ax.set_xlabel("Sequence position (bp; window midpoint)")
        ax.set_ylabel("GC content (%)")
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        output_paths.append(output_path)
        print(f"GC profile saved to: {output_path}")

    return output_paths


def plot_gc_profile(csv_path, output_path):
    """Backward-compatible function that plots the first sequence only."""
    csv_path = Path(csv_path)
    output_path = Path(output_path)

    # Use the multi-sequence plotting implementation, then retain only
    # the first plot at the requested legacy output path.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_dir = output_path.parent / "_gc_profile_temp"
    try:
        paths = plot_gc_profiles(csv_path, temporary_dir)

        first_plot = paths[0]
        output_path.write_bytes(first_plot.read_bytes())
        print(f"GC profile saved to: {output_path}")
    finally:
        if temporary_dir.exists():
            for path in temporary_dir.iterdir():
                if path.is_file():
                    path.unlink()
            temporary_dir.rmdir()


def main():
    parser = argparse.ArgumentParser(
        description="Plot sliding-window GC profiles for FASTA sequences."
    )
    parser.add_argument(
        "--input",
        default="output/gc_windows.csv",
        help="Path to the GC-window CSV",
    )
    parser.add_argument(
        "--output-dir",
        default="output/gc_profiles",
        help="Directory for the generated plots",
    )
    args = parser.parse_args()

    try:
        paths = plot_gc_profiles(args.input, args.output_dir)
        print(f"\nGenerated {len(paths)} GC profile(s).")
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

