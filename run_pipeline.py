

"""Run the complete SeqToolkit analysis workflow."""

import argparse
import csv
from pathlib import Path

from seqtoolkit import (
    parse_fasta,
    calculate_gc_windows,
)
from quality_control import quality_check
from plot_gc import plot_gc_profile


def main():
    parser = argparse.ArgumentParser(
        description="Run the SeqToolkit DNA analysis pipeline."
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--step-size", type=int, default=5)
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()

    if args.window_size < 1 or args.step_size < 1:
        parser.error("Window and step sizes must be positive.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        records = parse_fasta(args.input)
        if not records:
            parser.error("No sequences found in the FASTA file.")

        # 1. Quality-control report
        qc_results = [
            quality_check(header, sequence)
            for header, sequence in records
        ]

        qc_path = output_dir / "qc_report.csv"
        with qc_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=qc_results[0].keys()
            )
            writer.writeheader()
            writer.writerows(qc_results)

        # 2. Sliding-window GC report
        gc_rows = []
        for header, sequence in records:
            windows = calculate_gc_windows(
                sequence,
                window_size=args.window_size,
                step_size=args.step_size,
            )

            for window in windows:
                gc_rows.append({"header": header, **window})

        gc_path = output_dir / "gc_windows.csv"

        if gc_rows:
            with gc_path.open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=gc_rows[0].keys()
                )
                writer.writeheader()
                writer.writerows(gc_rows)
        else:
            # Avoid accidentally plotting a previous run's report.
            if gc_path.exists():
                gc_path.unlink()

        # 3. GC plot
        if gc_rows:
            plot_path = output_dir / "gc_profile.png"
            plot_gc_profile(gc_path, plot_path)
        else:
            print(
                "No complete GC windows were available; "
                "no GC plot was generated."
            )

        # 4. Summary
        print(f"\nProcessed {len(records)} FASTA record(s).")
        print(f"Quality report: {qc_path}")
        if gc_rows:
            print(f"GC-window report: {gc_path}")
            print(f"GC plot: {output_dir / 'gc_profile.png'}")

        print("\nQuality-control summary:")
        for result in qc_results:
            print(
                f'{result["status"]:7} | {result["header"]} | '
                f'length={result["length"]} | '
                f'ambiguous={result["ambiguous_count"]}'
            )

        print("\nPipeline completed.")

    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
