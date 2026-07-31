"""Tests for audited CAFE5 retry preparation after non-finite initialization."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREPARE_PATH = ROOT / "scripts" / "phylogeny" / "prepare_cafe5_timetree_bundle.py"
RETRY_PATH = ROOT / "scripts" / "phylogeny" / "prepare_cafe5_retry_bundle.py"
WORKFLOW_TEST_PATH = ROOT / "tests" / "phylogeny" / "test_cafe5_timetree_workflow.py"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


PREPARE = load("prepare_for_cafe_retry", PREPARE_PATH)
RETRY = load("prepare_cafe5_retry_bundle", RETRY_PATH)
WORKFLOW = load("cafe_workflow_retry_fixture", WORKFLOW_TEST_PATH)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Cafe5RetryPreparationTests(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, Path]:
        root = root.resolve()
        helper = WORKFLOW.Cafe5TimeTreeWorkflowTests()
        counts, dated, tree, executable = helper.make_fixture(root)
        bundle = root / "bundle"
        original = os.sys.argv
        try:
            os.sys.argv = [
                str(PREPARE_PATH), "--counts-dir", str(counts),
                "--dated-validation", str(dated / "validation.json"),
                "--dated-tree", str(tree), "--cafe5", str(executable),
                "--output-dir", str(bundle),
            ]
            self.assertEqual(PREPARE.main(), 0)
        finally:
            os.sys.argv = original
        checksums = digest(bundle / "checksums.tsv")
        failed = root / "failed"
        stage = failed / "runs/base_poisson"
        (stage / "results").mkdir(parents=True)
        (stage / "console.stdout").write_text(
            "Score (-lnL): inf\nFailed to initialize any reasonable values\n", encoding="utf-8"
        )
        (stage / "console.stderr").write_text(
            "Families with largest size differentials:\nOG2: 99\n", encoding="utf-8"
        )
        manifest = json.loads((bundle / "run_manifest.json").read_text())
        state = {
            "workflow": "sequential_cafe5_timetree_secondary",
            "status": "ERROR",
            "active_model": "base_poisson",
            "completed": [],
            "bundle": {
                "path": str(bundle),
                "manifest_sha256": digest(bundle / "run_manifest.json"),
                "checksums_sha256": checksums,
            },
            "error": str(stage / "results/Base_report.cafe") + " missing",
            "finished_at_utc": "2026-07-20T00:00:00+00:00",
            "cafe5": manifest["cafe5"],
        }
        (failed / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
        return bundle, failed

    @staticmethod
    def invoke(bundle: Path, failed: Path, output: Path) -> int:
        original = os.sys.argv
        try:
            os.sys.argv = [
                str(RETRY_PATH), "--source-bundle", str(bundle),
                "--failed-run", str(failed), "--output-dir", str(output),
            ]
            return RETRY.main()
        finally:
            os.sys.argv = original

    def test_retry_excludes_all_families_at_reported_difference_floor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, failed = self.fixture(root)
            output = root / "retry"
            self.assertEqual(self.invoke(bundle, failed, output), 0)
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["primary_family_count"], 1)
            self.assertEqual(manifest["runtime_outlier_count"], 1)
            self.assertEqual(
                manifest["retry_after_nonfinite_initialization"]["maximum_difference_exclusive"], 99
            )
            self.assertIn("OG2", (output / "cafe5_runtime_outliers.tsv").read_text())
            self.assertNotIn("OG2", (output / "cafe5_primary_lt100.tsv").read_text())

    def test_finite_failure_score_is_rejected_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, failed = self.fixture(root)
            (failed / "runs/base_poisson/console.stdout").write_text(
                "Score (-lnL): 123\nFailed to initialize any reasonable values\n", encoding="utf-8"
            )
            output = root / "retry"
            self.assertEqual(self.invoke(bundle, failed, output), 2)
            self.assertFalse(output.exists())

    def test_warning_difference_mismatch_is_rejected_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, failed = self.fixture(root)
            (failed / "runs/base_poisson/console.stderr").write_text("OG2: 98\n", encoding="utf-8")
            output = root / "retry"
            self.assertEqual(self.invoke(bundle, failed, output), 2)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
