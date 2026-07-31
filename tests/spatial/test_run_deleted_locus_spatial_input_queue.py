"""Tests for the detached deleted-locus spatial-input relay."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "run_deleted_locus_spatial_input_queue.py"


class DeletedLocusQueueTests(unittest.TestCase):
    def test_runs_once_after_both_prerequisites_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            search_state = root / "search.json"
            relabel_state = root / "relabel.json"
            search_state.write_text(
                json.dumps(
                    {
                        "workflow": "callable_aware_translated_search_queue",
                        "status": "PASS",
                        "completed": [{"unit": "u1"}],
                    }
                )
                + "\n"
            )
            relabel_state.write_text(
                json.dumps(
                    {
                        "workflow": "sequential_hy4a_similarity_relabelling_queue",
                        "status": "PASS",
                        "completed": [{"unit": "u1"}],
                    }
                )
                + "\n"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text("placeholder\n")
            worker = root / "worker.py"
            worker.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--manifest'); "
                "p.add_argument('--data-root'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "(o/'run_manifest.json').write_text(json.dumps({"
                "'status':'PASS','workflow':'callable_positive_deleted_expected_locus_spatial_inputs',"
                "'unit_count':1,'positive_deleted_count':2})+'\\n')\n"
            )
            output = root / "output"
            queue = root / "queue"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--search-queue-state",
                    str(search_state),
                    "--relabel-queue-state",
                    str(relabel_state),
                    "--manifest",
                    str(manifest),
                    "--data-root",
                    str(root),
                    "--output-dir",
                    str(output),
                    "--python",
                    sys.executable,
                    "--worker",
                    str(worker),
                    "--queue-root",
                    str(queue),
                    "--poll-seconds",
                    "1",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((queue / "state.json").read_text())
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["spatial_input_units"], 1)
            self.assertEqual(state["positive_deleted_count"], 2)


if __name__ == "__main__":
    unittest.main()
