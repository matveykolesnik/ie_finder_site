"""The Python driver: inputs, step logs, failures."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ie_finder import StepFailed, collect_inputs, step


class CollectInputsTest(unittest.TestCase):
    def test_directory_takes_fasta_suffixes_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name in ("a.fna", "b.fasta", "c.fa", "d.fsa", "notes.txt"):
                (tmp / name).write_text(">x\nACGT\n")
            self.assertEqual(sorted(collect_inputs(tmp)), ["a", "b", "c", "d"])

    def test_shared_sample_name_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "strain.fna").write_text(">x\nACGT\n")
            (tmp / "strain.fasta").write_text(">x\nACGT\n")
            with self.assertRaisesRegex(ValueError, "share the sample name strain"):
                collect_inputs(tmp)

    def test_sample_name_drops_the_last_suffix_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "T.oshimai.v2.fasta"
            path.write_text(">x\nACGT\n")
            self.assertEqual(list(collect_inputs(path)), ["T.oshimai.v2"])

    def test_empty_or_missing_input_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.fasta"
            empty.write_text("")
            with self.assertRaisesRegex(ValueError, "empty"):
                collect_inputs(empty)
            with self.assertRaisesRegex(ValueError, "not found"):
                collect_inputs(Path(tmp) / "missing.fasta")


class StepLogTest(unittest.TestCase):
    def test_python_and_child_process_output_land_in_the_log(self):
        def work():
            print("from python")
            subprocess.run(["sh", "-c", "echo from child >&2"], check=True)

        with tempfile.TemporaryDirectory() as tmp:
            step("S", Path(tmp), "demo", work)
            self.assertEqual((Path(tmp) / "demo.log").read_text(), "from python\nfrom child\n")

    def test_failure_names_the_step_and_logs_the_traceback(self):
        def work():
            raise RuntimeError("no tRNAs today")

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(StepFailed) as ctx:
                step("S", Path(tmp), "predict_trna", work)
            self.assertEqual((ctx.exception.sample, ctx.exception.step), ("S", "predict_trna"))
            self.assertIn("RuntimeError: no tRNAs today", ctx.exception.log.read_text())


if __name__ == "__main__":
    unittest.main()
