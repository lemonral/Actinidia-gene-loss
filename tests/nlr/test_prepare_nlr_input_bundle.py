"""Tests for atomic gzip/plain NLR input preparation."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "nlr" / "prepare_nlr_input_bundle.py"


class PrepareNlrInputBundleTest(unittest.TestCase):
    def write_fixture(self, root: Path) -> Path:
        (root / "reference.fa").write_text(">r1\nATGC\n", encoding="utf-8")
        with gzip.open(root / "target.fa.gz", "wb") as handle:
            handle.write(b">Chr01\nACGT\n>Chr02\nNNNN\n")
        manifest = root / "sources.tsv"
        manifest.write_text(
            "sample_id\tspecies\tploidy\tanalysis_role\tinput_scope\tsource_fasta\t"
            "output_basename\texpected_fasta_records\n"
            "reference\tClematoclethra scandens\tn/a\treference_callable\t"
            "reference_transcript_cds\treference.fa\treference.fa\t1\n"
            "target\tActinidia rufa\t2x\ttarget_repertoire\twhole_genome\t"
            "target.fa.gz\ttarget.fa\t2\n",
            encoding="utf-8",
        )
        return manifest

    def run_script(self, root: Path, manifest: Path, output_name: str = "bundle"):
        output = root / output_name
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--manifest",
                str(manifest),
                "--input-root",
                str(root),
                "--output-dir",
                str(output),
                "--expected-targets",
                "1",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        return completed, output

    def test_plain_and_gzip_inputs_are_materialized_and_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.write_fixture(root)
            completed, output = self.run_script(root, manifest)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual((output / "target.fa").read_bytes(), b">Chr01\nACGT\n>Chr02\nNNNN\n")
            with (output / "nlr_annotator_inputs.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["expected_fasta_records"] for row in rows], ["1", "2"])
            report = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["target_count"], 1)
            self.assertEqual({row["compression"] for row in report["sources"]}, {"plain", "gzip"})

    def test_bad_expected_record_count_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.write_fixture(root)
            text = manifest.read_text(encoding="utf-8").replace("target.fa\t2", "target.fa\t3")
            manifest.write_text(text, encoding="utf-8")
            completed, output = self.run_script(root, manifest)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("expected 3", completed.stderr)
            self.assertFalse(output.exists())

    def test_existing_output_and_duplicate_sample_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.write_fixture(root)
            output = root / "bundle"
            output.mkdir()
            completed, _ = self.run_script(root, manifest)
            self.assertEqual(completed.returncode, 2)
            manifest.write_text(
                manifest.read_text(encoding="utf-8").replace("target\tActinidia", "reference\tActinidia"),
                encoding="utf-8",
            )
            completed, duplicate_output = self.run_script(root, manifest, "duplicate")
            self.assertEqual(completed.returncode, 2)
            self.assertIn("duplicate sample_id", completed.stderr)
            self.assertFalse(duplicate_output.exists())


if __name__ == "__main__":
    unittest.main()
