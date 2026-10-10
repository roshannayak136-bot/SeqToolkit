
import csv
import tempfile
import unittest
from pathlib import Path

from seqtoolkit import (
    analyze_sequence,
    parse_fasta,
    write_csv,
)


class TestSequenceAnalysis(unittest.TestCase):

    def test_sequence_length_and_gc(self):
        result = analyze_sequence("test", "ATGC")
        self.assertEqual(result["length"], 4)
        self.assertEqual(result["GC_percent"], 50.0)

    def test_nucleotide_counts(self):
        result = analyze_sequence("test", "AACCGGTT")
        self.assertEqual(result["A_count"], 2)
        self.assertEqual(result["C_count"], 2)
        self.assertEqual(result["G_count"], 2)
        self.assertEqual(result["T_count"], 2)

    def test_reverse_complement(self):
        result = analyze_sequence("test", "ATGC")
        self.assertEqual(result["reverse_complement"], "GCAT")

    def test_ambiguous_bases(self):
        result = analyze_sequence("test", "ATGCRN")
        self.assertEqual(result["ambiguous_count"], 2)

    def test_invalid_symbol_is_rejected(self):
        with self.assertRaises(ValueError):
            analyze_sequence("test", "ATGX")

    def test_empty_sequence_is_rejected(self):
        with self.assertRaises(ValueError):
            analyze_sequence("test", "")

    def test_multiline_fasta(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta = Path(temp_dir) / "test.fasta"
            fasta.write_text(
                ">record1\nATG\nCA\n>record2\nTTTT\n",
                encoding="utf-8",
            )

            records = parse_fasta(fasta)

            self.assertEqual(
                records,
                [("record1", "ATGCA"), ("record2", "TTTT")],
            )

    def test_csv_export(self):
        result = analyze_sequence("test", "ATGC")

        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "summary.csv"
            write_csv([result], output)

            with output.open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["length"], "4")


    def test_gc_windows(self):
        from seqtoolkit import calculate_gc_windows

        results = calculate_gc_windows(
            "GGGGCCCCAAAA", window_size=4, step_size=4
        )

        self.assertEqual(
            results,
            [
                {"start": 1, "end": 4, "GC_percent": 100.0},
                {"start": 5, "end": 8, "GC_percent": 100.0},
                {"start": 9, "end": 12, "GC_percent": 0.0},
            ],
        )

    def test_gc_windows_reject_invalid_symbols(self):
        from seqtoolkit import calculate_gc_windows

        with self.assertRaises(ValueError):
            calculate_gc_windows("ATGX", window_size=2, step_size=1)

    def test_gc_windows_reject_zero_window(self):
        from seqtoolkit import calculate_gc_windows

        with self.assertRaises(ValueError):
            calculate_gc_windows("ATGC", window_size=0, step_size=1)

    def test_gc_windows_ignore_incomplete_tail(self):
        from seqtoolkit import calculate_gc_windows

        results = calculate_gc_windows(
            "ATGCG", window_size=4, step_size=4
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]["start"], 1)
        self.assertEqual(results[0]["end"], 4)

if __name__ == "__main__":
    unittest.main()
