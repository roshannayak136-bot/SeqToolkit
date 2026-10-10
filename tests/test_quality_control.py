
import unittest

from quality_control import quality_check


class TestQualityControl(unittest.TestCase):

    def test_unambiguous_sequence(self):
        result = quality_check("test", "ATGC")
        self.assertEqual(result["length"], 4)
        self.assertEqual(result["ambiguous_count"], 0)
        self.assertEqual(result["GC_percent_unambiguous"], 50.0)
        self.assertEqual(result["status"], "WARNING")

    def test_ambiguous_bases(self):
        result = quality_check("test", "ATGCNN")
        self.assertEqual(result["ambiguous_count"], 2)
        self.assertEqual(result["ambiguous_percent"], 33.33)
        self.assertEqual(result["GC_percent_unambiguous"], 50.0)
        self.assertIn("ambiguous", result["warnings"])

    def test_minimum_length_warning(self):
        result = quality_check("test", "ATGC", min_length=5)
        self.assertIn("shorter than minimum", result["warnings"])

    def test_maximum_length_warning(self):
        result = quality_check("test", "ATGCGC", max_length=5)
        self.assertIn("longer than maximum", result["warnings"])

    def test_invalid_symbol(self):
        with self.assertRaises(ValueError):
            quality_check("test", "ATGX")

    def test_empty_sequence(self):
        with self.assertRaises(ValueError):
            quality_check("test", "")

    def test_only_ambiguous_bases(self):
        result = quality_check("test", "NNNN")
        self.assertEqual(result["ambiguous_count"], 4)
        self.assertIsNone(result["GC_percent_unambiguous"])


if __name__ == "__main__":
    unittest.main()
