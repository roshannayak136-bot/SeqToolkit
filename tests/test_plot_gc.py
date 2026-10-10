
import csv
import tempfile
import unittest
from pathlib import Path

from plot_gc import plot_gc_profiles


class TestPlotGCProfiles(unittest.TestCase):

    def write_csv(self, path, rows, fieldnames=None):
        if fieldnames is None:
            fieldnames = ["header", "start", "end", "GC_percent"]

        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def test_creates_separate_plots_for_multiple_sequences(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            csv_path = root / "gc.csv"
            output_dir = root / "plots"

            self.write_csv(csv_path, [
                {
                    "header": "sequence_001",
                    "start": 0,
                    "end": 10,
                    "GC_percent": 50.0,
                },
                {
                    "header": "sequence_001",
                    "start": 5,
                    "end": 15,
                    "GC_percent": 60.0,
                },
                {
                    "header": "sequence_002",
                    "start": 0,
                    "end": 10,
                    "GC_percent": 40.0,
                },
            ])

            paths = plot_gc_profiles(csv_path, output_dir)

            self.assertEqual(len(paths), 2)
            self.assertTrue(all(path.exists() for path in paths))
            self.assertTrue(all(path.stat().st_size > 0 for path in paths))

    def test_empty_csv_raises_error(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            csv_path = root / "empty.csv"
            self.write_csv(csv_path, [])

            with self.assertRaisesRegex(ValueError, "No GC-window data"):
                plot_gc_profiles(csv_path, root / "plots")

    def test_missing_columns_raises_error(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            csv_path = root / "invalid.csv"

            self.write_csv(
                csv_path,
                [{"header": "seq1", "start": 0}],
                fieldnames=["header", "start"],
            )

            with self.assertRaisesRegex(ValueError, "missing required columns"):
                plot_gc_profiles(csv_path, root / "plots")

    def test_invalid_gc_percentage_raises_error(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            csv_path = root / "invalid_gc.csv"

            self.write_csv(csv_path, [{
                "header": "seq1",
                "start": 0,
                "end": 10,
                "GC_percent": 120,
            }])

            with self.assertRaisesRegex(ValueError, "Out-of-range"):
                plot_gc_profiles(csv_path, root / "plots")

    def test_invalid_window_coordinates_raise_error(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            csv_path = root / "invalid_window.csv"

            self.write_csv(csv_path, [{
                "header": "seq1",
                "start": 10,
                "end": 5,
                "GC_percent": 50,
            }])

            with self.assertRaisesRegex(ValueError, "Out-of-range"):
                plot_gc_profiles(csv_path, root / "plots")


if __name__ == "__main__":
    unittest.main()
