"""Regression tests for the historical shared-loss reproducer."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY_ROOT / "scripts" / "gene_loss" / "summarize_shared_loss.py"


class SummarizeSharedLossTest(unittest.TestCase):
    def test_summary_is_path_free_and_explicitly_historical(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            loss_dir = root / "legacy_loss_lists"
            loss_dir.mkdir()
            metadata = root / "samples.tsv"
            metadata.write_text(
                "target_haplotype\tspecies\tinclude_manuscript\n"
                "sample_a\tActinidia alpha\ttrue\n"
                "sample_b\tActinidia beta\ttrue\n",
                encoding="utf-8",
            )
            (loss_dir / "sample_a_decayed_genes.txt").write_text("gene_1\n", encoding="utf-8")
            (loss_dir / "sample_a_deleted_genes.txt").write_text("gene_2\n", encoding="utf-8")
            (loss_dir / "sample_b_decayed_genes.txt").write_text("gene_1\n", encoding="utf-8")
            (loss_dir / "sample_b_deleted_genes.txt").write_text("gene_3\n", encoding="utf-8")
            output = root / "output"

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--loss-dir",
                    str(loss_dir),
                    "--sample-metadata",
                    str(metadata),
                    "--cohort-name",
                    "historical_test",
                    "--run-id",
                    "test_v1",
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            summary_text = (output / "cohort_summary.json").read_text(encoding="utf-8")
            summary = json.loads(summary_text)
            self.assertEqual(summary["analysis_role"], "historical_manuscript_reproduction_only")
            self.assertEqual(summary["source_loss_directory_basename"], "legacy_loss_lists")
            self.assertNotIn("source_loss_directory", summary)
            self.assertNotIn(str(root), summary_text)
            self.assertEqual(summary["shared_loss_gene_count"], 1)


if __name__ == "__main__":
    unittest.main()
