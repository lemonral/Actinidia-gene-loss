"""Regression tests for biological-species gene-loss aggregation."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY_ROOT / "scripts" / "gene_loss" / "aggregate_species_loss.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AggregateSpeciesLossTest(unittest.TestCase):
    def write_metadata(self, root: Path) -> tuple[Path, list[str]]:
        haplotypes = ["Actinidia_eriantha_HAP1", "Actinidia_eriantha_HAP2"]
        subgenomes = [f"Actinidia_deliciosa_{suffix}" for suffix in "ABCDEF"]
        path = root / "species_units.tsv"
        lines = ["assembly_unit_id\tbiological_species\taggregation_rule\tinclude"]
        lines.extend(
            f"{unit}\tActinidia eriantha\tall_units_positive\ttrue" for unit in haplotypes
        )
        lines.extend(
            f"{unit}\tActinidia deliciosa\tany_unit_positive\ttrue" for unit in subgenomes
        )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path, haplotypes + subgenomes

    def run_script(self, root: Path, matrix_rows: list[tuple[str, str, str, str]]) -> subprocess.CompletedProcess[str]:
        metadata, _ = self.write_metadata(root)
        matrix = root / "unit_calls.tsv"
        with matrix.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["reference_gene_id", "assembly_unit_id", "classification", "callable"])
            writer.writerows(matrix_rows)
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--unit-call-matrix",
                str(matrix),
                "--unit-metadata",
                str(metadata),
                "--output-dir",
                str(root / "output"),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_haplotypes_and_subgenomes_are_aggregated_by_species(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            haplotypes = units[:2]
            subgenomes = units[2:]
            classifications: dict[str, dict[str, tuple[str, str]]] = {
                "gene_shared": {unit: ("pseudogenized", "true") for unit in units},
                "gene_partial": {
                    **{haplotypes[0]: ("deleted", "true"), haplotypes[1]: ("retained", "true")},
                    **{
                        unit: (("deleted", "true") if index == 0 else ("retained", "true"))
                        for index, unit in enumerate(subgenomes)
                    },
                },
                "gene_eriantha_only": {
                    **{unit: ("pseudogenized", "true") for unit in haplotypes},
                    **{unit: ("not_called_loss", "true") for unit in subgenomes},
                },
                "gene_uncertain": {
                    haplotypes[0]: ("retained", "true"),
                    haplotypes[1]: ("retained", "false"),
                    **{unit: ("retained", "true") for unit in subgenomes},
                },
            }
            rows = [
                (gene, unit, *classifications[gene][unit])
                for gene in classifications
                for unit in units
            ]
            completed = self.run_script(root, rows)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = root / "output"

            matrix = {
                (row["reference_gene_id"], row["biological_species"]): row
                for row in read_tsv(output / "species_gene_matrix.tsv")
            }
            self.assertEqual(
                matrix[("gene_shared", "Actinidia eriantha")]["species_gene_status"],
                "positive_complete",
            )
            self.assertEqual(
                matrix[("gene_shared", "Actinidia deliciosa")]["species_gene_status"],
                "positive_complete",
            )
            self.assertEqual(
                matrix[("gene_partial", "Actinidia eriantha")]["species_gene_status"],
                "positive_partial",
            )
            self.assertEqual(
                matrix[("gene_partial", "Actinidia eriantha")]["species_positive_by_rule"],
                "false",
            )
            self.assertEqual(
                matrix[("gene_partial", "Actinidia deliciosa")]["positive_unit_count"], "1"
            )
            self.assertEqual(
                matrix[("gene_partial", "Actinidia deliciosa")]["species_positive_by_rule"],
                "true",
            )
            self.assertEqual(
                matrix[("gene_uncertain", "Actinidia eriantha")]["species_gene_status"],
                "uncertain",
            )
            self.assertEqual(
                matrix[("gene_uncertain", "Actinidia eriantha")]["uncertain_unit_count"], "1"
            )

            shared = read_tsv(output / "shared_positive_complete_genes.tsv")
            self.assertEqual([row["reference_gene_id"] for row in shared], ["gene_shared"])
            non_shared = read_tsv(output / "non_shared_positive_calls.tsv")
            partial = [row for row in non_shared if row["reference_gene_id"] == "gene_partial"]
            self.assertEqual(len(partial), 1)
            self.assertEqual(partial[0]["biological_species"], "Actinidia deliciosa")
            self.assertEqual(partial[0]["confident_lineage_restricted_species_loss"], "false")
            complete = [
                row for row in non_shared if row["reference_gene_id"] == "gene_eriantha_only"
            ]
            self.assertEqual(len(complete), 1)
            self.assertEqual(complete[0]["confident_lineage_restricted_species_loss"], "false")

            unit_long = read_tsv(output / "unit_calls_long.tsv")
            self.assertEqual(len(unit_long), len(classifications) * len(units))
            summary = json.loads((output / "species_loss_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["schema_version"], "2.0")
            self.assertEqual(summary["status"], "PASS")
            self.assertEqual(summary["assembly_unit_count"], 8)
            self.assertEqual(summary["biological_species_count"], 2)
            self.assertEqual(summary["shared_positive_complete_gene_count"], 1)
            self.assertTrue(all(summary["checks"].values()))
            self.assertEqual(
                [row["biological_species"] for row in summary["species_aggregation"]],
                ["Actinidia deliciosa", "Actinidia eriantha"],
            )
            serialized = json.dumps(summary)
            self.assertNotIn(str(root), serialized)
            self.assertEqual({item["basename"] for item in summary["inputs"]}, {"unit_calls.tsv", "species_units.tsv"})

    def test_zero_positive_cohort_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            rows = [("gene_zero", unit, "retained", "true") for unit in units]
            completed = self.run_script(root, rows)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = root / "output"
            self.assertEqual(read_tsv(output / "shared_positive_complete_genes.tsv"), [])
            self.assertEqual(read_tsv(output / "non_shared_positive_calls.tsv"), [])
            prevalence = read_tsv(output / "species_prevalence.tsv")
            self.assertEqual(len(prevalence), 1)
            self.assertEqual(prevalence[0]["positive_by_rule_species_count"], "0")
            summary = json.loads((output / "species_loss_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["shared_positive_complete_gene_count"], 0)
            self.assertEqual(summary["non_shared_positive_call_count"], 0)

    def test_duplicate_unit_gene_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            rows = [("gene_1", unit, "retained", "true") for unit in units]
            rows.append(("gene_1", units[0], "retained", "true"))
            completed = self.run_script(root, rows)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("duplicate unit-gene row", completed.stderr)
            self.assertFalse((root / "output").exists())

    def test_missing_unit_gene_grid_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            rows = [("gene_1", unit, "retained", "true") for unit in units]
            rows.extend(("gene_2", unit, "retained", "true") for unit in units[:-1])
            completed = self.run_script(root, rows)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("incomplete selected unit-by-gene grid", completed.stderr)
            self.assertFalse((root / "output").exists())

    def test_species_rule_must_be_explicit_and_consistent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata = root / "species_units.tsv"
            metadata.write_text(
                "assembly_unit_id\tbiological_species\taggregation_rule\tinclude\n"
                "Actinidia_eriantha_HAP1\tActinidia eriantha\tall_units_positive\ttrue\n"
                "Actinidia_eriantha_HAP2\tActinidia eriantha\tany_unit_positive\ttrue\n",
                encoding="utf-8",
            )
            matrix = root / "unit_calls.tsv"
            matrix.write_text(
                "reference_gene_id\tassembly_unit_id\tclassification\tcallable\n"
                "gene_1\tActinidia_eriantha_HAP1\tretained\ttrue\n"
                "gene_1\tActinidia_eriantha_HAP2\tretained\ttrue\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--unit-call-matrix",
                    str(matrix),
                    "--unit-metadata",
                    str(metadata),
                    "--output-dir",
                    str(root / "output"),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("inconsistent aggregation rules", completed.stderr)

    def test_positive_noncallable_is_rejected_and_not_called_is_uncertain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            rows = [("gene_1", unit, "retained", "true") for unit in units]
            rows[0] = ("gene_1", units[0], "deleted", "false")
            completed = self.run_script(root, rows)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("positive classification 'deleted' requires callable=true", completed.stderr)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            _, units = self.write_metadata(root)
            rows = [("gene_1", unit, "not_called_loss", "true") for unit in units]
            completed = self.run_script(root, rows)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            matrix = read_tsv(root / "output" / "species_gene_matrix.tsv")
            self.assertTrue(all(row["species_gene_status"] == "uncertain" for row in matrix))


if __name__ == "__main__":
    unittest.main()
