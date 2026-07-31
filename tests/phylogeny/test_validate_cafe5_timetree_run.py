"""Tests for fail-closed validation of TimeTree-bound CAFE5 output."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREPARE_PATH = ROOT / "scripts" / "phylogeny" / "prepare_cafe5_timetree_bundle.py"
VALIDATE_PATH = ROOT / "scripts" / "phylogeny" / "validate_cafe5_timetree_run.py"
WORKFLOW_TEST_PATH = ROOT / "tests" / "phylogeny" / "test_cafe5_timetree_workflow.py"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


PREPARE = load("prepare_for_cafe_validation", PREPARE_PATH)
VALIDATE = load("validate_cafe5_timetree_run", VALIDATE_PATH)
WORKFLOW = load("cafe_workflow_test_fixture", WORKFLOW_TEST_PATH)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Cafe5ValidationTests(unittest.TestCase):
    def prepare_bundle(self, root: Path) -> Path:
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
        return bundle

    @staticmethod
    def write_model(results: Path, prefix: str, family_ids: list[str]) -> None:
        results.mkdir(parents=True)
        nodes = ["A<0>", "B<1>", "C<2>", "<3>", "<4>"]
        header = ["FamilyID", *nodes]
        for suffix, rows in (
            ("count.tab", [[family, "1", "2", "3", "2", "2"] for family in family_ids]),
            ("change.tab", [[family, "0", "1", "1", "0", "0"] for family in family_ids]),
        ):
            with (results / f"{prefix}_{suffix}").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)
        (results / f"{prefix}_family_results.txt").write_text(
            "#FamilyID\tpvalue\tSignificant at 0.05\n" + "".join(
                f"{family}\t{0.01 if index == 0 else 0.5}\t{'y' if index == 0 else 'n'}\n"
                for index, family in enumerate(family_ids)
            ), encoding="utf-8",
        )
        (results / f"{prefix}_clade_results.txt").write_text(
            "#Taxon_ID\tIncrease\tDecrease\nA<0>\t1\t0\n", encoding="utf-8"
        )
        if prefix == "Base":
            likelihoods = "#FamilyID\tLikelihood of Family\n" + "".join(
                f"{family}\t0.5\n" for family in family_ids
            )
        else:
            likelihoods = "#FamilyID\tGamma Cat Mean\tLkhd of Category\tLkhd of Family\tPosterior Probability\tSignificant\n" + "".join(
                f"{family}\t1\t0.5\t0.5\t1\tN/S\n" for family in family_ids
            )
        (results / f"{prefix}_family_likelihoods.txt").write_text(likelihoods, encoding="utf-8")
        (results / f"{prefix}_branch_probabilities.tab").write_text(
            "\t".join(header) + "\n" + "".join(
                f"{family}\t0.5\t0.5\t0.5\t0.5\tN/A\n" for family in family_ids
            ), encoding="utf-8",
        )
        (results / f"{prefix}_asr.tre").write_text(
            "#nexus\nBEGIN TREES;\n" + "".join(
                f" TREE {family} = ((A_1:50,B_2:50)_2:50,C_3:100)_2;\n" for family in family_ids
            ) + "END;\n", encoding="utf-8",
        )
        (results / f"{prefix}_results.txt").write_text(
            f"Model {prefix} Result: 123.5\nLambda: 0.002\n", encoding="utf-8"
        )
        (results / f"{prefix}_report.cafe").write_text(
            "Tree:((A:50,B:50):50,C:100):0\nLambda:\t0.002\nLambda tree:\t((1,1)1,1)1\n",
            encoding="utf-8",
        )

    def make_valid_run(self, root: Path, bundle: Path) -> Path:
        root = root.resolve()
        bundle = bundle.resolve()
        manifest = json.loads((bundle / "run_manifest.json").read_text())
        family_ids = ["OG1", "OG2"]
        run = root / "run"
        (run / "runs").mkdir(parents=True)
        completed = []
        for model_id, prefix, arguments in (
            ("base_poisson", "Base", ["--poisson"]),
            ("gamma3_poisson", "Gamma", ["--poisson", "--n_gamma_cats", "3"]),
        ):
            stage = run / "runs" / model_id
            results = stage / "results"
            self.write_model(results, prefix, family_ids)
            (stage / "console.stdout").write_text("complete\n", encoding="utf-8")
            (stage / "console.stderr").write_text("", encoding="utf-8")
            inventory = sorted(path for path in results.iterdir() if path.is_file())
            completed.append({
                "model_id": model_id,
                "returncode": 0,
                "started_at_utc": "2026-07-20T00:00:00+00:00",
                "finished_at_utc": "2026-07-20T00:01:00+00:00",
                "command": [
                    manifest["cafe5"]["path"], "--infile", str(bundle / "cafe5_primary_lt100.tsv"),
                    "--tree", str(bundle / "dated_tree.mean_ma.tre"), "--cores", "1",
                    *arguments, "--output_prefix", str(results),
                ],
                "stdout": {"bytes": (stage / "console.stdout").stat().st_size, "sha256": digest(stage / "console.stdout")},
                "stderr": {"bytes": 0, "sha256": digest(stage / "console.stderr")},
                "result_files": [
                    {"relative_path": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}
                    for path in inventory
                ],
            })
        state = {
            "schema_version": 1,
            "workflow": "sequential_cafe5_timetree_secondary",
            "status": "PASS_RUN_COMPLETE",
            "calibration_claim": VALIDATE.CALIBRATION_CLAIM,
            "bundle": {
                "path": str(bundle),
                "manifest_sha256": digest(bundle / "run_manifest.json"),
                "checksums_sha256": digest(bundle / "checksums.tsv"),
            },
            "cafe5": manifest["cafe5"],
            "cores": 1,
            "large_family_count_excluded_from_rate_estimation": manifest["large_family_count"],
            "runtime_outlier_count_excluded_from_rate_estimation": manifest.get(
                "runtime_outlier_count", 0
            ),
            "completed": completed,
        }
        (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
        return run

    @staticmethod
    def invoke(
        bundle: Path, run: Path, output: Path, *, accept_gamma_failure: bool = False,
    ) -> int:
        original = os.sys.argv
        try:
            os.sys.argv = [
                str(VALIDATE_PATH), "--bundle", str(bundle), "--run-dir", str(run),
                "--output-dir", str(output),
            ]
            if accept_gamma_failure:
                os.sys.argv.append("--accept-gamma-initialization-failure")
            return VALIDATE.main()
        finally:
            os.sys.argv = original

    def test_valid_run_is_atomic_and_summarized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            output = root / "validation"
            self.assertEqual(self.invoke(bundle, run, output), 0)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "PASS_CAFE5_VALIDATED")
            self.assertEqual(validation["primary_family_count"], 2)
            self.assertIn("base_poisson\tOG1\t0.01", (output / "significant_families.tsv").read_text())
            self.assertTrue((output / "checksums.tsv").is_file())

    def test_running_state_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            state = json.loads((run / "state.json").read_text())
            state["status"] = "running"
            (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
            output = root / "validation"
            self.assertEqual(self.invoke(bundle, run, output), 2)
            self.assertFalse(output.exists())

    def test_result_tampering_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            (run / "runs/base_poisson/results/Base_results.txt").write_text("tampered\n")
            output = root / "validation"
            self.assertEqual(self.invoke(bundle, run, output), 2)
            self.assertFalse(output.exists())

    def test_family_closure_failure_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            path = run / "runs/base_poisson/results/Base_family_results.txt"
            path.write_text("#FamilyID\tpvalue\tSignificant at 0.05\nOG1\t0.01\ty\n", encoding="utf-8")
            state = json.loads((run / "state.json").read_text())
            binding = next(
                item for item in state["completed"][0]["result_files"]
                if item["relative_path"] == path.name
            )
            binding.update(bytes=path.stat().st_size, sha256=digest(path))
            (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
            output = root / "validation"
            self.assertEqual(self.invoke(bundle, run, output), 2)
            self.assertFalse(output.exists())

    def test_exact_gamma_initialization_failure_validates_base_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            state = json.loads((run / "state.json").read_text())
            gamma_stage = run / "runs/gamma3_poisson"
            for path in (gamma_stage / "results").iterdir():
                path.unlink()
            (gamma_stage / "results").rmdir()
            (gamma_stage / "console.stdout").write_text(
                "Starting Search for Initial Parameter Values\n"
                "Failed to initialize any reasonable values\n",
                encoding="utf-8",
            )
            (gamma_stage / "console.stderr").write_text(
                "Families with largest size differentials:\nOG1: 10\n\n"
                "You may want to try removing the top few families\n",
                encoding="utf-8",
            )
            state["status"] = "ERROR"
            state["completed"] = state["completed"][:1]
            state["error"] = (
                "missing or symlink file: "
                + str(gamma_stage / "results" / "Gamma_report.cafe")
            )
            (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
            output = root / "validation"
            self.assertEqual(
                self.invoke(bundle, run, output, accept_gamma_failure=True), 0
            )
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(
                validation["status"], "PASS_CAFE5_BASE_VALIDATED_GAMMA_UNAVAILABLE"
            )
            self.assertEqual(validation["models"][0]["model_id"], "base_poisson")
            self.assertEqual(
                validation["unavailable_sensitivity"]["status"],
                "UNAVAILABLE_INITIALIZATION_FAILURE",
            )

    def test_gamma_failure_mode_rejects_incomplete_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = self.prepare_bundle(root)
            run = self.make_valid_run(root, bundle)
            state = json.loads((run / "state.json").read_text())
            gamma_stage = run / "runs/gamma3_poisson"
            for path in (gamma_stage / "results").iterdir():
                path.unlink()
            (gamma_stage / "results").rmdir()
            (gamma_stage / "console.stdout").write_text(
                "Failed to initialize any reasonable values\n", encoding="utf-8"
            )
            (gamma_stage / "console.stderr").write_text("truncated\n", encoding="utf-8")
            state["status"] = "ERROR"
            state["completed"] = state["completed"][:1]
            state["error"] = (
                "missing or symlink file: "
                + str(gamma_stage / "results" / "Gamma_report.cafe")
            )
            (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
            output = root / "validation"
            self.assertEqual(
                self.invoke(bundle, run, output, accept_gamma_failure=True), 2
            )
            self.assertFalse(output.exists())

    def test_declared_root_filter_is_bound_as_analyzed_subset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "Base_count.tab"
            table.write_text(
                "FamilyID\tA<0>\nOG1\t1\nOG3\t1\n", encoding="utf-8"
            )
            console = root / "console.stdout"
            console.write_text(
                "Filtering families not present at the root from: 3 to 2\n",
                encoding="utf-8",
            )
            self.assertEqual(
                VALIDATE.read_analyzed_family_ids(table, ["OG1", "OG2", "OG3"], console),
                ["OG1", "OG3"],
            )

    def test_root_filter_rejects_unbound_subset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "Base_count.tab"
            table.write_text("FamilyID\tA<0>\nOG1\t1\n", encoding="utf-8")
            console = root / "console.stdout"
            console.write_text("no filtering declaration\n", encoding="utf-8")
            with self.assertRaises(VALIDATE.ValidationError):
                VALIDATE.read_analyzed_family_ids(table, ["OG1", "OG2"], console)


if __name__ == "__main__":
    unittest.main()
