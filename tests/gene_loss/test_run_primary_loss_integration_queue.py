"""Tests for the gated primary-matrix merge and aggregation relay."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "run_primary_loss_integration_queue.py"


class PrimaryLossIntegrationQueueTests(unittest.TestCase):
    def test_merges_then_aggregates_after_upstream_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            upstream = root / "upstream.json"
            upstream.write_text(
                json.dumps(
                    {
                        "workflow": "postsearch_complete_loss_matrix_and_aggregation_queue",
                        "status": "PASS",
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            sources = root / "sources.tsv"
            metadata = root / "metadata.tsv"
            sources.write_text("placeholder\n", encoding="utf-8")
            metadata.write_text("placeholder\n", encoding="utf-8")
            merge_worker = root / "merge.py"
            merge_worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--sources'); p.add_argument('--data-root'); "
                "p.add_argument('--expected-total-units'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'complete_unit_loss_matrix.tsv').write_text('reference_gene_id\\n'); "
                "(o/'run_manifest.json').write_text(json.dumps({"
                "'status':'PASS','workflow':'merged_primary_complete_loss_matrix',"
                "'assembly_unit_count':2,'reference_gene_count':3,'matrix_rows':6,'expected_matrix_rows':6})+'\\n')\n",
                encoding="utf-8",
            )
            aggregate_worker = root / "aggregate.py"
            aggregate_worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--unit-call-matrix'); "
                "p.add_argument('--unit-metadata'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'species_loss_summary.json').write_text(json.dumps({"
                "'status':'PASS','assembly_unit_count':2,'reference_gene_count':3,"
                "'biological_species_count':2,'shared_positive_complete_gene_count':1})+'\\n')\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--upstream-state", str(upstream),
                    "--sources", str(sources), "--data-root", str(root),
                    "--expected-total-units", "2", "--merged-output", str(root / "merged"),
                    "--aggregation-metadata", str(metadata),
                    "--aggregation-output", str(root / "aggregation"),
                    "--python", sys.executable, "--merge-worker", str(merge_worker),
                    "--aggregation-worker", str(aggregate_worker),
                    "--queue-root", str(root / "queue"), "--poll-seconds", "1",
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((root / "queue" / "state.json").read_text())
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["assembly_unit_count"], 2)
            self.assertEqual(state["biological_species_count"], 2)


if __name__ == "__main__":
    unittest.main()
