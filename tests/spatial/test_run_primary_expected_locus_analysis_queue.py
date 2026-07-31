"""Tests for the gated primary expected-locus spatial analysis relay."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "run_primary_expected_locus_analysis_queue.py"


class PrimaryExpectedLocusQueueTests(unittest.TestCase):
    def test_merges_and_analyzes_after_both_inputs_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            new_state = root / "new.json"
            new_state.write_text(json.dumps({"workflow": "deleted_locus_spatial_input_queue", "status": "PASS"}))
            legacy = root / "legacy"
            legacy.mkdir()
            (legacy / "run_manifest.json").write_text(
                json.dumps({"workflow": "legacy_conservative_deleted_expected_locus_spatial_inputs", "status": "PASS"})
            )
            sources = root / "sources.tsv"
            sources.write_text("placeholder\n")
            merge = root / "merge.py"
            merge.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--sources'); p.add_argument('--data-root'); "
                "p.add_argument('--expected-total-units'); p.add_argument('--output-dir'); a=p.parse_args()\n"
                "o=pathlib.Path(a.output_dir); o.mkdir(); "
                "[(o/n).write_text('placeholder\\n') for n in ['positive_deleted_calls.tsv','expected_deleted_locus_coordinates.tsv','assembly_manifest.tsv']]; "
                "(o/'run_manifest.json').write_text(json.dumps({'status':'PASS','workflow':'merged_primary_expected_deleted_locus_spatial_inputs','unit_count':2,'positive_deleted_count':3})+'\\n')\n"
            )
            analysis = root / "analysis.py"
            analysis.write_text(
                "import argparse,json,pathlib\n"
                "p=argparse.ArgumentParser(); [p.add_argument(x) for x in ['--positive-calls','--feature-coordinates','--assembly-manifest','--output-dir','--analysis-label','--positive-classes','--number-of-bins','--call-unit-column','--call-gene-column','--call-classification-column','--coordinate-unit-column','--coordinate-gene-column','--coordinate-chromosome-column','--coordinate-start-column','--coordinate-end-column','--coordinate-classification-column']]; p.add_argument('--legacy-reproduction',action='store_true'); a=p.parse_args(); "
                "o=pathlib.Path(a.output_dir); o.mkdir(); (o/'run_summary.json').write_text(json.dumps({'reconciliation':{'assembly_unit_count':2,'positive_call_count':3,'emitted_position_count':3},'centromere_policy':'none'})+'\\n')\n"
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--new-spatial-state", str(new_state),
                    "--legacy-input-dir", str(legacy), "--sources", str(sources),
                    "--data-root", str(root), "--expected-total-units", "2",
                    "--merged-input-dir", str(root / "merged"),
                    "--analysis-output-dir", str(root / "result"), "--python", sys.executable,
                    "--merge-worker", str(merge), "--analysis-worker", str(analysis),
                    "--queue-root", str(root / "queue"), "--poll-seconds", "1",
                ], text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((root / "queue" / "state.json").read_text())
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["positive_deleted_count"], 3)


if __name__ == "__main__":
    unittest.main()
