"""Tests for the gated callable-copy-opportunity summary relay."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "run_copy_opportunity_summary_queue.py"


class CopyOpportunitySummaryQueueTests(unittest.TestCase):
    def test_runs_after_primary_integration_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            upstream = root / "upstream.json"
            upstream.write_text(
                json.dumps(
                    {
                        "workflow": "primary_complete_loss_integration_queue", "status": "PASS",
                        "assembly_unit_count": 2, "biological_species_count": 1,
                        "reference_gene_count": 3,
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            aggregation = root / "aggregation"
            aggregation.mkdir()
            for name in ("unit_calls_long.tsv", "species_gene_matrix.tsv", "shared_positive_complete_genes.tsv"):
                (aggregation / name).write_text("placeholder\n", encoding="utf-8")
            worker = root / "worker.py"
            worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--unit-calls'); p.add_argument('--species-matrix'); "
                "p.add_argument('--shared-genes'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'run_manifest.json').write_text(json.dumps({"
                "'status':'PASS','workflow':'callable_copy_opportunity_and_loss_mode_summary',"
                "'assembly_unit_count':2,'biological_species_count':1,'reference_gene_count':3,"
                "'shared_positive_complete_gene_count':1})+'\\n')\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--upstream-state", str(upstream),
                    "--aggregation-dir", str(aggregation), "--output-dir", str(root / "out"),
                    "--python", sys.executable, "--worker", str(worker),
                    "--queue-root", str(root / "queue"), "--poll-seconds", "1",
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((root / "queue" / "state.json").read_text())
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["shared_positive_complete_gene_count"], 1)


if __name__ == "__main__":
    unittest.main()
