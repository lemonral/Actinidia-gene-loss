"""Tests for TimeTree-bound CAFE5 preparation and sequential execution."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREPARE_SCRIPT = ROOT / "scripts" / "phylogeny" / "prepare_cafe5_timetree_bundle.py"
RUN_SCRIPT = ROOT / "scripts" / "phylogeny" / "run_cafe5_timetree_bundle.py"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


PREPARE = load("prepare_cafe5_timetree_bundle", PREPARE_SCRIPT)
RUN = load("run_cafe5_timetree_bundle", RUN_SCRIPT)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_checksums(directory: Path, names: list[str]) -> None:
    (directory / "checksums.tsv").write_text(
        "file\tsha256\n" + "".join(f"{name}\t{digest(directory / name)}\n" for name in names),
        encoding="utf-8",
    )


class Cafe5TimeTreeWorkflowTests(unittest.TestCase):
    def make_fixture(self, root: Path, *, tree_text: str = "((A:50,B:50):50,C:100);"):
        counts = root / "counts"
        counts.mkdir()
        matrix = counts / "cafe5_family_counts.tsv"
        matrix.write_text(
            "Desc\tFamily ID\tA\tB\tC\n"
            "NA\tOG1\t1\t2\t3\n"
            "NA\tOG2\t99\t0\t1\n"
            "NA\tOG3\t100\t1\t0\n",
            encoding="utf-8",
        )
        summary = {
            "schema_version": 1,
            "workflow": "topology_bound_cafe5_family_count_preparation",
            "status": "PASS_INPUT_PREPARATION_ONLY",
            "matrix": {"basename": matrix.name, "bytes": matrix.stat().st_size, "sha256": digest(matrix)},
            "terminal_order": ["A", "B", "C"],
            "family_count": 3,
            "families_with_any_terminal_count_at_least_100": 1,
        }
        (counts / "preparation.json").write_text(json.dumps(summary) + "\n", encoding="utf-8")
        write_checksums(counts, [matrix.name, "preparation.json"])

        dated = root / "dated"
        dated.mkdir()
        tree = dated / "dated_tree.mean_ma.tre"
        tree.write_text(tree_text + "\n", encoding="ascii")
        validation = {
            "status": "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC",
            "calibration_claim": PREPARE.CALIBRATION_CLAIM,
            "dated_tree": {"root_age_ma": 100.0, "maximum_root_to_tip_deviation_ma": 0.0},
        }
        (dated / "validation.json").write_text(json.dumps(validation) + "\n", encoding="utf-8")
        write_checksums(dated, [tree.name, "validation.json"])

        install = root / "CAFE5_5.1.0"
        (install / "bin").mkdir(parents=True)
        executable = install / "bin" / "cafe5"
        executable.write_text(
            """#!/usr/bin/env python3
import pathlib, sys
if '--help' in sys.argv:
    print('Usage: cafe5 [options]\\n--infile')
    raise SystemExit(0)
prefix = pathlib.Path(sys.argv[sys.argv.index('--output_prefix') + 1])
prefix.mkdir()
model = 'Gamma' if '--n_gamma_cats' in sys.argv else 'Base'
for suffix in ('report.cafe','results.txt','family_likelihoods.txt','asr.tre','count.tab','change.tab','family_results.txt','clade_results.txt','branch_probabilities.tab'):
    (prefix / f'{model}_{suffix}').write_text(f'{model} output\\n')
""",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        files = {
            "CHANGELOG.md": "## [1.1.0]\n",
            "README.md": "CAFE 5\n",
            "LICENSE": "license\n",
            "main.cpp": "int main() {}\n",
            "Makefile": "CXX=g++\n",
            "config.log": "PACKAGE_VERSION='1.1'\n",
        }
        for name, text in files.items():
            (install / name).write_text(text, encoding="utf-8")
        return counts, dated, tree, executable

    def prepare(self, counts: Path, dated: Path, tree: Path, executable: Path, output: Path) -> int:
        original = os.sys.argv
        try:
            os.sys.argv = [
                str(PREPARE_SCRIPT), "--counts-dir", str(counts),
                "--dated-validation", str(dated / "validation.json"),
                "--dated-tree", str(tree), "--cafe5", str(executable),
                "--output-dir", str(output),
            ]
            return PREPARE.main()
        finally:
            os.sys.argv = original

    def invoke_runner(self, bundle: Path, output: Path) -> int:
        original = os.sys.argv
        try:
            os.sys.argv = [str(RUN_SCRIPT), "--bundle", str(bundle), "--output-dir", str(output)]
            return RUN.main()
        finally:
            os.sys.argv = original

    def test_preparation_splits_large_families_and_binds_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts, dated, tree, executable = self.make_fixture(root)
            bundle = root / "bundle"
            self.assertEqual(self.prepare(counts, dated, tree, executable, bundle), 0)
            manifest = json.loads((bundle / "run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "PASS_PREPARED_CAFE5")
            self.assertEqual(manifest["primary_family_count"], 2)
            self.assertEqual(manifest["large_family_count"], 1)
            self.assertIn("OG3", (bundle / "cafe5_large_ge100.tsv").read_text())
            self.assertNotIn("OG3", (bundle / "cafe5_primary_lt100.tsv").read_text())

    def test_non_ultrametric_tree_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts, dated, tree, executable = self.make_fixture(
                root, tree_text="((A:50,B:49):50,C:100);"
            )
            bundle = root / "bundle"
            self.assertEqual(self.prepare(counts, dated, tree, executable, bundle), 2)
            self.assertFalse(bundle.exists())

    def test_runner_completes_base_and_gamma_sequentially(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts, dated, tree, executable = self.make_fixture(root)
            bundle = root / "bundle"
            self.assertEqual(self.prepare(counts, dated, tree, executable, bundle), 0)
            output = root / "run"
            self.assertEqual(self.invoke_runner(bundle, output), 0)
            state = json.loads((output / "state.json").read_text())
            self.assertEqual(state["status"], "PASS_RUN_COMPLETE")
            self.assertEqual([row["model_id"] for row in state["completed"]], [
                "base_poisson", "gamma3_poisson"
            ])
            self.assertTrue((output / "runs/base_poisson/results/Base_report.cafe").is_file())
            self.assertTrue((output / "runs/gamma3_poisson/results/Gamma_report.cafe").is_file())

    def test_runner_rejects_bundle_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts, dated, tree, executable = self.make_fixture(root)
            bundle = root / "bundle"
            self.assertEqual(self.prepare(counts, dated, tree, executable, bundle), 0)
            (bundle / "cafe5_primary_lt100.tsv").write_text("tampered\n", encoding="utf-8")
            output = root / "run"
            self.assertEqual(self.invoke_runner(bundle, output), 2)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
