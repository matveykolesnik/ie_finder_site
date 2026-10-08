"""Base search parameters with an override on top."""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from merge_search_params import merge_search_params

BASE = ROOT / "search_params.yaml"


class MergeSearchParamsTest(unittest.TestCase):
    def _write(self, text):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "strict.yaml"
        path.write_text(text)
        return path

    def test_override_keys_win_and_others_keep_base_values(self):
        merged = merge_search_params(BASE, self._write("attl_select_by: length\n"))
        self.assertEqual(merged.attl_select_by, "length")
        self.assertEqual(merged.attl_window_bp, 300000)

    def test_no_override_gives_the_base(self):
        self.assertEqual(merge_search_params(BASE, None).attl_select_by, "bitscore")

    def test_error_names_the_override_file(self):
        path = self._write("attl_windw_bp: 1\n")
        with self.assertRaisesRegex(ValueError, "strict.yaml: unknown search parameter"):
            merge_search_params(BASE, path)


if __name__ == "__main__":
    unittest.main()
