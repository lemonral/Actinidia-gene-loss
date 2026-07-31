"""Regression and fail-closed tests for biological-species PGLS."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY_ROOT / "scripts" / "downstream" / "species_pgls.py"
DEPENDENCIES_AVAILABLE = all(
    importlib.util.find_spec(name) is not None for name in ("numpy", "scipy", "Bio")
)

sys.path.insert(0, str(REPOSITORY_ROOT / "src"))
import geneloss_repro.pgls as pgls_module  # noqa: E402


def binding(path: Path) -> dict[str, str | int]:
    payload = path.read_bytes()
    return {
        "basename": path.name,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PGLSDependencyContractTest(unittest.TestCase):
    def test_scientific_dependencies_must_not_be_silently_skipped(self) -> None:
        self.assertTrue(
            DEPENDENCIES_AVAILABLE,
            "Run PGLS tests with an interpreter containing NumPy, SciPy, and Biopython; "
            "install the project with python -m pip install -e '.[statistics]'",
        )


@unittest.skipUnless(DEPENDENCIES_AVAILABLE, "species PGLS optional dependencies unavailable")
class SpeciesPGLSTest(unittest.TestCase):
    species = [f"Species {letter}" for letter in "ABCDEF"]

    def write_tree(self, root: Path, text: str | None = None) -> Path:
        path = root / "time_tree.nwk"
        path.write_text(
            text
            or "((('Species A':1,'Species B':1):1,('Species C':1,'Species D':1):1):1,"
            "('Species E':2,'Species F':2):1);\n",
            encoding="utf-8",
        )
        return path

    def rows(self) -> list[dict[str, str]]:
        counts = [50, 77, 110, 151, 205, 280]
        return [
            {
                "biological_species": species,
                "analysis_level": "biological_species",
                "loss_scope": "lineage_specific_nonshared",
                "lineage_specific_nonshared_positive_loss_count": str(count),
                "callable_denominator": "1000",
                "log2_ploidy": str(index),
            }
            for index, (species, count) in enumerate(zip(self.species, counts), start=1)
        ]

    def write_data(
        self,
        root: Path,
        rows: list[dict[str, str]] | None = None,
        extra_field: str | None = None,
    ) -> Path:
        path = root / "species.tsv"
        records = rows or self.rows()
        fields = list(records[0])
        if extra_field and extra_field not in fields:
            fields.append(extra_field)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
        return path

    def write_pass_reports(
        self,
        root: Path,
        *,
        data: Path,
        tree: Path,
        species: list[str] | None = None,
    ) -> dict[str, Path]:
        expected_species = species or self.species
        matrix = root / "species_gene_matrix.tsv"
        matrix.write_text("validated species matrix fixture\n", encoding="utf-8")
        shared = root / "shared_positive_complete_genes.tsv"
        shared.write_text("reference_gene_id\n", encoding="utf-8")
        ledger = root / "ploidy_ledger.tsv"
        ledger.write_text("validated ploidy ledger fixture\n", encoding="utf-8")
        species_loss = root / "species_loss_summary.json"
        species_loss_document = {
            "schema_version": "2.0",
            "status": "PASS",
            "definitions": {},
            "inputs": [],
            "include_column": "include_gene_loss",
            "assembly_unit_count": len(expected_species),
            "biological_species_count": len(expected_species),
            "reference_gene_count": 1000,
            "expected_unit_matrix_rows": 1000 * len(expected_species),
            "species_gene_matrix_rows": 1000 * len(expected_species),
            "species_gene_status_counts": {"not_positive": 1000 * len(expected_species)},
            "shared_positive_complete_gene_count": 0,
            "non_shared_positive_call_count": 0,
            "confident_lineage_restricted_gene_count": 0,
            "aggregation_rule_species_counts": {
                "all_units_positive": len(expected_species)
            },
            "species_aggregation": [
                {
                    "biological_species": species_name,
                    "aggregation_rule": "all_units_positive",
                    "assembly_unit_count": 1,
                    "assembly_units": [f"unit_{index}"],
                }
                for index, species_name in enumerate(expected_species, start=1)
            ],
            "checks": {
                "complete_selected_unit_gene_grid": True,
                "positive_classification_requires_callable": True,
                "not_called_loss_treated_as_uncertain": True,
                "species_status_counts_reconciled": True,
                "shared_set_reconciled": True,
                "output_checksums_reconciled": True,
            },
            "outputs": [],
        }
        species_loss.write_text(
            json.dumps(species_loss_document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        ploidy_report = root / "ploidy_ledger_pass.json"
        ploidy_report.write_text(
            json.dumps(
                {
                    "schema_version": "species_ploidy_ledger_pass_v1",
                    "workflow": "species_ploidy_ledger_validation",
                    "workflow_version": "1.0.0",
                    "status": "PASS",
                    "analysis_level": "biological_species",
                    "predictor": "log2_ploidy",
                    "ploidy_ledger": binding(ledger),
                    "biological_species": expected_species,
                    "checks": {
                        "one_row_per_biological_species": True,
                        "positive_integer_ploidy": True,
                        "log2_ploidy_recalculated_exactly": True,
                        "exact_biological_species_set": True,
                        "source_provenance_complete": True,
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        input_report = root / "pgls_input_pass.json"
        input_report.write_text(
            json.dumps(
                {
                    "schema_version": "species_pgls_input_pass_v1",
                    "workflow": "species_pgls_input_builder",
                    "workflow_version": "1.0.0",
                    "status": "PASS",
                    "analysis_level": "biological_species",
                    "loss_scope": "lineage_specific_nonshared",
                    "predictor": "log2_ploidy",
                    "input_data": binding(data),
                    "species_count": len(expected_species),
                    "biological_species": expected_species,
                    "aggregation_policy": pgls_module.AGGREGATION_POLICY,
                    "upstream_bindings": {
                        "species_loss_manifest": binding(species_loss),
                        "species_gene_matrix": binding(matrix),
                        "shared_positive_complete_gene_set": binding(shared),
                        "ploidy_ledger": binding(ledger),
                        "ploidy_ledger_pass_report": binding(ploidy_report),
                    },
                    "checks": {key: True for key in sorted(pgls_module.INPUT_PASS_CHECKS)},
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        dating_manifest = root / "dating_manifest.json"
        dating_manifest.write_text('{"status":"PASS"}\n', encoding="utf-8")
        tree_report = root / "time_tree_pass.json"
        tree_report.write_text(
            json.dumps(
                {
                    "schema_version": "species_time_tree_pass_v1",
                    "workflow": "species_time_tree_validation",
                    "workflow_version": "1.0.0",
                    "status": "PASS",
                    "analysis_level": "biological_species",
                    "tree": binding(tree),
                    "source_dating_manifest": binding(dating_manifest),
                    "biological_species": expected_species,
                    "root_semantics": "accepted_biological_species_mrca",
                    "branch_length_units": "million_years",
                    "checks": {key: True for key in sorted(pgls_module.TREE_PASS_CHECKS)},
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return {
            "input": input_report,
            "species_loss": species_loss,
            "ploidy": ploidy_report,
            "tree": tree_report,
        }

    def run_analysis(
        self,
        root: Path,
        *,
        data: Path | None = None,
        tree: Path | None = None,
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        data_path = data or self.write_data(root)
        tree_path = tree or self.write_tree(root)
        reports = self.write_pass_reports(root, data=data_path, tree=tree_path)
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--data",
                str(data_path),
                "--time-tree",
                str(tree_path),
                "--input-pass-report",
                str(reports["input"]),
                "--species-loss-manifest",
                str(reports["species_loss"]),
                "--ploidy-ledger-pass-report",
                str(reports["ploidy"]),
                "--time-tree-pass-report",
                str(reports["tree"]),
                "--predictor-column",
                "log2_ploidy",
                "--output-dir",
                str(root / "output"),
                *(extra_args or []),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_complete_bundle_primary_loo_and_named_sensitivity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            completed = self.run_analysis(
                root,
                extra_args=["--sensitivity", "without_f=Species F"],
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = root / "output"
            expected = {
                "analysis_data.tsv",
                "model_summary.tsv",
                "model_coefficients.tsv",
                "fitted_residuals.tsv",
                "leave_one_species_out.tsv",
                "named_exclusion_sensitivities.tsv",
                "publication_gate.tsv",
                "analysis_manifest.json",
                "checksums.sha256.tsv",
            }
            self.assertEqual({path.name for path in output.iterdir()}, expected)
            summaries = read_tsv(output / "model_summary.tsv")
            self.assertEqual([row["model_id"] for row in summaries], ["primary", "sensitivity:without_f"])
            self.assertTrue(all(0 <= float(row["lambda_ml"]) <= 1 for row in summaries))
            coefficients = read_tsv(output / "model_coefficients.tsv")
            predictor = next(
                row for row in coefficients if row["model_id"] == "primary" and row["term"] == "log2_ploidy"
            )
            self.assertGreater(float(predictor["estimate"]), 0)
            self.assertGreater(float(predictor["standard_error"]), 0)
            self.assertEqual(len(read_tsv(output / "fitted_residuals.tsv")), 6)
            self.assertEqual(len(read_tsv(output / "leave_one_species_out.tsv")), 6)
            sensitivities = read_tsv(output / "named_exclusion_sensitivities.tsv")
            self.assertEqual(sensitivities[0]["excluded_species"], "Species F")
            manifest_text = (output / "analysis_manifest.json").read_text(encoding="utf-8")
            manifest = json.loads(manifest_text)
            self.assertEqual(manifest["schema_version"], "species_pgls_v2")
            self.assertEqual(manifest["loss_scope"], "lineage_specific_nonshared")
            self.assertEqual(manifest["response_formula"], "log((L+0.5)/(D-L+0.5))")
            self.assertEqual(
                manifest["publication_gate"],
                "BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED",
            )
            self.assertNotIn(str(root), manifest_text)
            gate = read_tsv(output / "publication_gate.tsv")
            self.assertEqual(
                gate[0]["publication_status"],
                "BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED",
            )
            checksums = read_tsv(output / "checksums.sha256.tsv")
            self.assertEqual({row["relative_path"] for row in checksums}, expected - {"checksums.sha256.tsv"})

    def assert_fails_without_output(
        self,
        root: Path,
        expected_message: str,
        *,
        rows: list[dict[str, str]] | None = None,
        tree_text: str | None = None,
        extra_field: str | None = None,
        extra_args: list[str] | None = None,
    ) -> None:
        data = self.write_data(root, rows=rows, extra_field=extra_field)
        tree = self.write_tree(root, tree_text)
        completed = self.run_analysis(root, data=data, tree=tree, extra_args=extra_args)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(expected_message, completed.stderr)
        self.assertFalse((root / "output").exists())

    def test_duplicate_species_and_technical_unit_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            rows[-1]["biological_species"] = "Species E"
            self.assert_fails_without_output(root, "duplicate biological species", rows=rows)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            for row in rows:
                row["assembly_unit_id"] = row["biological_species"] + "_HAP1"
            self.assert_fails_without_output(
                root,
                "technical-unit columns are forbidden",
                rows=rows,
                extra_field="assembly_unit_id",
            )

    def test_scope_level_and_count_contracts_are_rejected(self) -> None:
        mutations = [
            ("loss_scope", "all_positive_losses", "must be exactly 'lineage_specific_nonshared'"),
            ("analysis_level", "haplotype", "must be exactly 'biological_species'"),
            (
                "lineage_specific_nonshared_positive_loss_count",
                "1001",
                "exceeds callable_denominator",
            ),
            (
                "lineage_specific_nonshared_positive_loss_count",
                "-1",
                "must be a non-negative integer",
            ),
            (
                "lineage_specific_nonshared_positive_loss_count",
                "1.5",
                "must be a non-negative integer",
            ),
            ("callable_denominator", "0", "must be positive"),
        ]
        for field, value, message in mutations:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary_directory:
                root = Path(temporary_directory)
                rows = self.rows()
                rows[0][field] = value
                self.assert_fails_without_output(root, message, rows=rows)

    def test_exact_tip_reconciliation_and_ultrametricity_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = "((('Species X':1,'Species B':1):1,('Species C':1,'Species D':1):1):1,('Species E':2,'Species F':2):1);"
            self.assert_fails_without_output(root, "exact tree/data tip reconciliation failed", tree_text=tree)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = "((('Species A':1.2,'Species B':1):1,('Species C':1,'Species D':1):1):1,('Species E':2,'Species F':2):1);"
            self.assert_fails_without_output(root, "time tree is not ultrametric", tree_text=tree)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = "((('Species A','Species B':1):1,('Species C':1,'Species D':1):1):1,('Species E':2,'Species F':2):1);"
            self.assert_fails_without_output(
                root, "every non-root branch needs a finite, non-negative length", tree_text=tree
            )

    def test_singular_covariance_and_insufficient_n_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            tree = "((('Species A':0,'Species B':0):2,('Species C':1,'Species D':1):1):1,('Species E':2,'Species F':2):1);"
            self.assert_fails_without_output(root, "Brownian covariance is singular", tree_text=tree)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()[:5]
            tree = "((('Species A':1,'Species B':1):1,('Species C':1,'Species D':1):1):1,'Species E':3);"
            self.assert_fails_without_output(root, "requires at least 6", rows=rows, tree_text=tree)

    def test_constant_predictor_and_invalid_named_exclusion_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            for row in rows:
                row["log2_ploidy"] = "1"
            self.assert_fails_without_output(root, "has no variation", rows=rows)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            for row in rows:
                row["lineage_specific_nonshared_positive_loss_count"] = "10"
            self.assert_fails_without_output(root, "response has no variation", rows=rows)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root,
                "species absent from data",
                extra_args=["--sensitivity", "bad=Species Z"],
            )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root,
                "requires at least 5 biological species",
                extra_args=["--sensitivity", "too_small=Species A,Species B"],
            )

    def test_exact_large_integer_parsing_and_noncanonical_numbers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            exact_denominator = "9007199254740993"
            for row in rows:
                row["callable_denominator"] = exact_denominator
            data = self.write_data(root, rows=rows)
            observations = pgls_module.read_species_data(
                data, predictor_column="log2_ploidy"
            )
            self.assertTrue(
                all(item.callable_denominator == 9007199254740993 for item in observations)
            )
        for invalid in ("1e3", "1_000", "+1000", "1.0"):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as temporary_directory:
                root = Path(temporary_directory)
                rows = self.rows()
                rows[0]["callable_denominator"] = invalid
                data = self.write_data(root, rows=rows)
                with self.assertRaisesRegex(pgls_module.SchemaError, "non-negative integer"):
                    pgls_module.read_species_data(data, predictor_column="log2_ploidy")
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            rows[0]["callable_denominator"] = str(2**63)
            data = self.write_data(root, rows=rows)
            with self.assertRaisesRegex(pgls_module.SchemaError, "signed-64-bit limit"):
                pgls_module.read_species_data(data, predictor_column="log2_ploidy")
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = self.rows()
            rows[0]["log2_ploidy"] = "1_0"
            data = self.write_data(root, rows=rows)
            with self.assertRaisesRegex(pgls_module.SchemaError, "canonical non-negative decimal"):
                pgls_module.read_species_data(data, predictor_column="log2_ploidy")

    def test_duplicate_header_and_ragged_row_are_rejected_without_traceback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            text = data.read_text(encoding="utf-8")
            header, rest = text.split("\n", 1)
            data.write_text(
                header.replace(
                    "callable_denominator",
                    "callable_denominator\tcallable_denominator",
                )
                + "\n"
                + rest.replace("\t1000\t", "\t1000\t1000\t"),
                encoding="utf-8",
            )
            completed = self.run_analysis(root, data=data)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("duplicate TSV header", completed.stderr)
            self.assertNotIn("Traceback", completed.stderr)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            lines = data.read_text(encoding="utf-8").splitlines()
            lines[1] += "\textra"
            data.write_text("\n".join(lines) + "\n", encoding="utf-8")
            completed = self.run_analysis(root, data=data)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("ragged TSV row", completed.stderr)
            self.assertNotIn("Traceback", completed.stderr)

    def test_input_mutation_after_snapshot_blocks_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            tree = self.write_tree(root)
            reports = self.write_pass_reports(root, data=data, tree=tree)
            original = pgls_module.build_brownian_covariance

            def mutate_source_then_build(*args, **kwargs):
                data.write_text(
                    data.read_text(encoding="utf-8").replace("\t50\t1000\t1\n", "\t51\t1000\t1\n"),
                    encoding="utf-8",
                )
                return original(*args, **kwargs)

            with mock.patch.object(
                pgls_module,
                "build_brownian_covariance",
                side_effect=mutate_source_then_build,
            ):
                with self.assertRaisesRegex(pgls_module.SchemaError, "input changed after validation"):
                    pgls_module.run_species_pgls(
                        data_path=data,
                        tree_path=tree,
                        input_pass_report_path=reports["input"],
                        species_loss_manifest_path=reports["species_loss"],
                        ploidy_ledger_pass_report_path=reports["ploidy"],
                        tree_pass_report_path=reports["tree"],
                        output_dir=root / "output",
                        predictor_column="log2_ploidy",
                    )
            self.assertFalse((root / "output").exists())

    def test_atomic_no_replace_preserves_racing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            tree = self.write_tree(root)
            reports = self.write_pass_reports(root, data=data, tree=tree)
            original = pgls_module._rename_directory_no_replace

            def inject_destination(source: Path, destination: Path) -> None:
                destination.mkdir()
                (destination / "sentinel.txt").write_text("keep\n", encoding="utf-8")
                original(source, destination)

            with mock.patch.object(
                pgls_module,
                "_rename_directory_no_replace",
                side_effect=inject_destination,
            ):
                with self.assertRaisesRegex(pgls_module.SchemaError, "refusing overwrite"):
                    pgls_module.run_species_pgls(
                        data_path=data,
                        tree_path=tree,
                        input_pass_report_path=reports["input"],
                        species_loss_manifest_path=reports["species_loss"],
                        ploidy_ledger_pass_report_path=reports["ploidy"],
                        tree_pass_report_path=reports["tree"],
                        output_dir=root / "output",
                        predictor_column="log2_ploidy",
                    )
            self.assertEqual((root / "output" / "sentinel.txt").read_text(), "keep\n")

    def test_relative_ultrametric_tolerance_polytomy_and_star_fail_closed(self) -> None:
        tiny = (
            "((('Species A':0.0000001,'Species B':0.0000002):0,"
            "('Species C':0.0000001,'Species D':0.0000002):0):0,"
            "('Species E':0.0000001,'Species F':0.0000002):0);"
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root, "time tree is not ultrametric", tree_text=tiny
            )
        polytomy = (
            "(('Species A':1,'Species B':1,'Species C':1):1,"
            "(('Species D':1,'Species E':1):1,'Species F':2):0);"
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root, "must be strictly bifurcating", tree_text=polytomy
            )
        star = ",".join(f"'{species}':1" for species in self.species)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root, "must be strictly bifurcating", tree_text=f"({star});"
            )
        ill_conditioned = (
            "(((\'Species A\':1e-14,\'Species B\':1e-14):1,"
            "(\'Species C\':0.5,\'Species D\':0.5):0.5):1,"
            "(\'Species E\':1,\'Species F\':1):1);"
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self.assert_fails_without_output(
                root, "condition number", tree_text=ill_conditioned
            )

    def test_flat_pagel_lambda_profile_is_rejected(self) -> None:
        import numpy as np

        observations = [
            pgls_module.SpeciesObservation(
                species,
                count,
                1000,
                float(index),
                float(__import__("math").log((count + 0.5) / (1000 - count + 0.5))),
            )
            for index, (species, count) in enumerate(
                zip(self.species, [50, 77, 110, 151, 205, 280]), start=1
            )
        ]
        with self.assertRaisesRegex(pgls_module.SchemaError, "profile likelihood is flat"):
            pgls_module.fit_pgls(observations, np.eye(6), model_id="star")

    def test_upstream_pass_and_all_units_positive_contracts_are_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            tree = self.write_tree(root)
            reports = self.write_pass_reports(root, data=data, tree=tree)
            manifest = json.loads(reports["species_loss"].read_text(encoding="utf-8"))
            manifest["species_aggregation"][0]["aggregation_rule"] = "any_unit_positive"
            reports["species_loss"].write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(pgls_module.SchemaError, "requires all_units_positive"):
                pgls_module.run_species_pgls(
                    data_path=data,
                    tree_path=tree,
                    input_pass_report_path=reports["input"],
                    species_loss_manifest_path=reports["species_loss"],
                    ploidy_ledger_pass_report_path=reports["ploidy"],
                    tree_pass_report_path=reports["tree"],
                    output_dir=root / "output",
                    predictor_column="log2_ploidy",
                )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data = self.write_data(root)
            tree = self.write_tree(root)
            reports = self.write_pass_reports(root, data=data, tree=tree)
            report = json.loads(reports["ploidy"].read_text(encoding="utf-8"))
            report["status"] = "BLOCKED"
            reports["ploidy"].write_text(
                json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(pgls_module.SchemaError, "invalid status"):
                pgls_module.run_species_pgls(
                    data_path=data,
                    tree_path=tree,
                    input_pass_report_path=reports["input"],
                    species_loss_manifest_path=reports["species_loss"],
                    ploidy_ledger_pass_report_path=reports["ploidy"],
                    tree_pass_report_path=reports["tree"],
                    output_dir=root / "output",
                    predictor_column="log2_ploidy",
                )


if __name__ == "__main__":
    unittest.main()
