"""Tests for the detached primary NLR downstream relay."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
QUEUE = ROOT / "scripts" / "nlr" / "run_primary_nlr_summary_queue.py"
DIGESTS = ("1" * 64, "2" * 64, "3" * 64)


def write_script(path: Path, source: str) -> None:
    path.write_text(textwrap.dedent(source), encoding="utf-8")


class RunPrimaryNlrSummaryQueueTest(unittest.TestCase):
    def fixture(self, root: Path, *, published: bool = True) -> dict[str, Path]:
        upstream = root / "nlr"
        if published:
            upstream.mkdir()
        bundle = root / "bundle"
        bundle.mkdir()
        for name in ("units.tsv", "loss.tsv", "shared.tsv"):
            (root / name).write_text("header\nvalue\n", encoding="utf-8")
        (root / "controller.pid").write_text("999999\n", encoding="utf-8")
        adapter = root / "adapter.py"
        write_script(
            adapter,
            """
            import json, pathlib, sys
            out = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])
            out.mkdir()
            for name in ('assembly_units.tsv', 'repertoire_counts.tsv',
                         'positive_reference_nlr_loss_calls.tsv',
                         'callable_reference_nlr_denominators.tsv'):
                (out / name).write_text('header\\nvalue\\n')
            (out / 'validation.json').write_text(json.dumps({
                'status': 'PASS_PRIMARY_NLR_SUMMARY_INPUTS',
                'assembly_unit_count': 23, 'reference_gene_count': 35547,
                'shared_positive_complete_gene_count': 68,
                'reference_nlr_gene_count': 300,
                'nonshared_reference_nlr_gene_count': 295,
            }))
            """,
        )
        summarizer = root / "summarizer.py"
        write_script(
            summarizer,
            """
            import json, pathlib, sys
            out = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])
            out.mkdir()
            (out / 'nlr_unit_summary.tsv').write_text('header\\nvalue\\n')
            (out / 'validation.json').write_text(json.dumps({
                'status': 'pass', 'analysis_cohort': 'primary_23_units_nonshared_v1',
                'cohort_role': 'primary', 'assembly_unit_count': 23,
                'denominator_input_mode': 'catalog',
                'positive_reference_nlr_loss_call_count': 50,
                'callable_reference_nlr_denominator_sum': 5000,
            }))
            """,
        )
        figure = root / "figure.py"
        write_script(
            figure,
            """
            import hashlib, json, pathlib, sys
            out = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])
            name = sys.argv[sys.argv.index('--basename') + 1]
            out.mkdir()
            payloads = {
                name + '.png': b'png', name + '.pdf': b'pdf',
                name + '.plot_data.tsv': b'data', name + '.caption.txt': b'caption',
                name + '.validation.json': json.dumps({
                    'status': 'pass', 'assembly_unit_count': 23
                }).encode(),
            }
            outputs = []
            for basename, content in payloads.items():
                path = out / basename
                path.write_bytes(content)
                outputs.append({
                    'basename': basename, 'bytes': len(content),
                    'sha256': hashlib.sha256(content).hexdigest(),
                })
            (out / (name + '.manifest.json')).write_text(json.dumps({
                'bundle_basename': name, 'outputs': outputs,
            }))
            """,
        )
        return {
            "upstream": upstream, "bundle": bundle, "metadata": root / "units.tsv",
            "loss": root / "loss.tsv", "shared": root / "shared.tsv",
            "pid": root / "controller.pid", "adapter": adapter,
            "summarizer": summarizer, "figure": figure,
        }

    def command(self, root: Path, values: dict[str, Path]) -> list[str]:
        return [
            sys.executable, str(QUEUE),
            "--upstream-controller-pid-file", str(values["pid"]),
            "--nlr-root", str(values["upstream"]), "--input-bundle", str(values["bundle"]),
            "--unit-metadata", str(values["metadata"]), "--loss-matrix", str(values["loss"]),
            "--shared-positive-genes", str(values["shared"]),
            "--input-output-dir", str(root / "inputs"), "--summary-output-dir",
            str(root / "summary"), "--figure-output-dir", str(root / "figure_output"),
            "--queue-root", str(root / "queue"), "--python", sys.executable,
            "--adapter", str(values["adapter"]), "--summarizer", str(values["summarizer"]),
            "--figure-worker", str(values["figure"]), "--source-root", str(ROOT / "src"),
            "--expected-jar-sha256", DIGESTS[0], "--expected-motifs-sha256", DIGESTS[1],
            "--expected-store-sha256", DIGESTS[2], "--poll-seconds", "1",
        ]

    def test_completed_upstream_runs_three_workers_sequentially_and_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            values = self.fixture(root)
            completed = subprocess.run(
                self.command(root, values), cwd=ROOT, capture_output=True, text=True, check=False
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            state = json.loads((root / "queue" / "state.json").read_text(encoding="utf-8"))
            self.assertEqual(state["status"], "PASS")
            self.assertEqual(state["reference_nlr_gene_count"], 300)
            self.assertEqual(state["positive_reference_nlr_loss_call_count"], 50)
            self.assertEqual(state["figure_output_checksum_count"], 5)
            self.assertTrue((root / "figure_output" / "primary_nlr_nonshared.png").is_file())

    def test_dead_upstream_before_atomic_publication_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            values = self.fixture(root, published=False)
            completed = subprocess.run(
                self.command(root, values), cwd=ROOT, capture_output=True, text=True, check=False
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("exited before atomic publication", completed.stderr)
            state = json.loads((root / "queue" / "state.json").read_text(encoding="utf-8"))
            self.assertEqual(state["status"], "ERROR")
            self.assertFalse((root / "inputs").exists())


if __name__ == "__main__":
    unittest.main()
