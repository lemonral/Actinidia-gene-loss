"""Tests for the post-search complete-matrix aggregation relay."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "run_postsearch_loss_aggregation_queue.py"


class PostsearchAggregationQueueTests(unittest.TestCase):
    def test_builds_matrix_then_aggregation_after_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            search = root / "search.json"
            search.write_text(
                json.dumps(
                    {
                        "workflow": "callable_aware_translated_search_queue",
                        "status": "PASS",
                    }
                )
                + "\n"
            )
            manifest = root / "matrix.tsv"
            metadata = root / "metadata.tsv"
            manifest.write_text("placeholder\n")
            metadata.write_text("placeholder\n")
            matrix_worker = root / "matrix_worker.py"
            matrix_worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--manifest'); p.add_argument('--data-root'); "
                "p.add_argument('--reference-cds'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'complete_unit_loss_matrix.tsv').write_text('reference_gene_id\\n'); "
                "(o/'run_manifest.json').write_text(json.dumps({"
                "'status':'PASS','workflow':'complete_callable_aware_new_unit_loss_matrix',"
                "'assembly_unit_count':2,'reference_gene_count':3,'matrix_rows':6,'expected_matrix_rows':6})+'\\n')\n"
            )
            aggregate_worker = root / "aggregate_worker.py"
            aggregate_worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--unit-call-matrix'); "
                "p.add_argument('--unit-metadata'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'species_loss_summary.json').write_text(json.dumps({"
                "'status':'PASS','assembly_unit_count':2,'reference_gene_count':3,"
                "'biological_species_count':1,'shared_positive_complete_gene_count':1})+'\\n')\n"
            )
            matrix_output = root / "matrix_output"
            aggregate_output = root / "aggregate_output"
            queue = root / "queue"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--search-queue-state", str(search),
                    "--matrix-manifest", str(manifest),
                    "--reference-cds", "reference.fa",
                    "--data-root", str(root),
                    "--matrix-output", str(matrix_output),
                    "--aggregation-metadata", str(metadata),
                    "--aggregation-output", str(aggregate_output),
                    "--python", sys.executable,
                    "--matrix-worker", str(matrix_worker),
                    "--aggregation-worker", str(aggregate_worker),
                    "--queue-root", str(queue),
                    "--poll-seconds", "1",
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((queue / "state.json").read_text())
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["matrix_rows"], 6)
            self.assertEqual(state["biological_species_count"], 1)
            self.assertEqual(state["shared_positive_complete_gene_count"], 1)


if __name__ == "__main__":
    unittest.main()
