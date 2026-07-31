"""Tests for the detached species-PGLS relay."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "run_species_pgls_queue.py"
SPEC = importlib.util.spec_from_file_location("run_species_pgls_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class SpeciesPGLSQueueTests(unittest.TestCase):
    def test_passes_exact_fake_workers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            species = [f"Species {letter}" for letter in "ABCDEF"]
            upstream = root / "upstream.json"
            upstream.write_text(json.dumps({
                "workflow": "primary_complete_loss_integration_queue", "status": "PASS"
            }) + "\n")
            loss = root / "loss"; loss.mkdir()
            (loss / "species_loss_summary.json").write_text("{}\n")
            ploidy = root / "ploidy.tsv"; ploidy.write_text("fixture\n")
            tree_dir = root / "tree"; tree_dir.mkdir()
            (tree_dir / "species_time_tree.nwk").write_text("(A:1,B:1):0;\n")
            (tree_dir / "species_time_tree_pass.json").write_text(json.dumps({
                "status": "PASS", "workflow": "species_time_tree_validation",
                "biological_species": species,
            }) + "\n")
            builder = root / "builder.py"
            builder.write_text(
                "import argparse,hashlib,json,pathlib\n"
                "p=argparse.ArgumentParser(); p.add_argument('--species-loss-dir'); p.add_argument('--ploidy-ledger'); p.add_argument('--output-dir'); a=p.parse_args(); o=pathlib.Path(a.output_dir); o.mkdir()\n"
                f"(o/'pgls_input_pass.json').write_text(json.dumps({{'status':'PASS','workflow':'species_pgls_input_builder','biological_species':{species!r}}})+'\\n')\n"
                "(o/'pgls_input.tsv').write_text('fixture\\n'); (o/'ploidy_ledger_pass.json').write_text('{}\\n')\n"
                "fs=sorted(x for x in o.iterdir()); (o/'checksums.sha256.tsv').write_text('relative_path\\tsha256\\n'+''.join(f'{x.name}\\t{hashlib.sha256(x.read_bytes()).hexdigest()}\\n' for x in fs))\n"
            )
            worker = root / "worker.py"
            worker.write_text(
                "import argparse,hashlib,json,pathlib\n"
                "p=argparse.ArgumentParser(); [p.add_argument(x) for x in ['--data','--time-tree','--input-pass-report','--species-loss-manifest','--ploidy-ledger-pass-report','--time-tree-pass-report','--predictor-column','--sensitivity','--output-dir']]; a=p.parse_args(); o=pathlib.Path(a.output_dir); o.mkdir()\n"
                f"m={{'status':'COMPLETE_EXPLORATORY_BLOCKED_FOR_PUBLICATION','publication_gate':'BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED','species_count':{len(species)},'named_exclusion_sensitivities':{{'without_rufa':['Actinidia rufa']}}}}\n"
                "(o/'analysis_manifest.json').write_text(json.dumps(m)+'\\n'); x=o/'analysis_manifest.json'; (o/'checksums.sha256.tsv').write_text('relative_path\\tsha256\\n'+f'{x.name}\\t{hashlib.sha256(x.read_bytes()).hexdigest()}\\n')\n"
            )
            old = sys.argv
            input_dir, output_dir, queue = root / "input", root / "output", root / "queue"
            try:
                sys.argv = [
                    str(SCRIPT), "--upstream-state", str(upstream), "--species-loss-dir", str(loss),
                    "--ploidy-ledger", str(ploidy), "--time-tree-dir", str(tree_dir),
                    "--input-dir", str(input_dir), "--output-dir", str(output_dir),
                    "--python", sys.executable, "--input-builder", str(builder),
                    "--pgls-worker", str(worker), "--queue-root", str(queue), "--poll-seconds", "1",
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                sys.argv = old
            state = json.loads((queue / "state.json").read_text())
            self.assertEqual(state["status"], "PASS_EXPLORATORY_PGLS")


if __name__ == "__main__":
    unittest.main()
