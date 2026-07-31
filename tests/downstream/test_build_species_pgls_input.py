"""Tests for the checksum-bound species-PGLS input builder."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AGGREGATE = ROOT / "scripts" / "gene_loss" / "aggregate_species_loss.py"
BUILD = ROOT / "scripts" / "downstream" / "build_species_pgls_input.py"
SPECIES = [f"Actinidia species_{index}" for index in range(1, 7)]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BuildSpeciesPGLSInputTest(unittest.TestCase):
    def make_aggregation(self, root: Path, *, any_unit: bool = False) -> Path:
        metadata = root / "metadata.tsv"
        with metadata.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly_unit_id", "biological_species", "aggregation_rule", "include"])
            for index, species in enumerate(SPECIES, start=1):
                rule = "any_unit_positive" if any_unit and index == 1 else "all_units_positive"
                writer.writerow([f"unit_{index}", species, rule, "true"])
        matrix = root / "unit_matrix.tsv"
        with matrix.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["reference_gene_id", "assembly_unit_id", "classification", "callable"])
            for species_index in range(1, 7):
                writer.writerow(["g_shared", f"unit_{species_index}", "deleted", "true"])
            for gene_index in range(1, 7):
                for species_index in range(1, 7):
                    classification = "deleted" if gene_index <= species_index else "retained"
                    writer.writerow(
                        [f"g_{gene_index}", f"unit_{species_index}", classification, "true"]
                    )
            for species_index in range(1, 7):
                writer.writerow(
                    ["g_uncertain", f"unit_{species_index}", "not_called_loss", "true"]
                )
        output = root / "species_loss"
        completed = subprocess.run(
            [
                sys.executable,
                str(AGGREGATE),
                "--unit-call-matrix",
                str(matrix),
                "--unit-metadata",
                str(metadata),
                "--output-dir",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return output

    def write_ploidy(self, root: Path, species: list[str] | None = None) -> Path:
        path = root / "ploidy.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["biological_species", "ploidy", "ploidy_source", "source_reference"])
            for index, name in enumerate(species or SPECIES, start=1):
                writer.writerow([name, [2, 2, 4, 4, 6, 6][index - 1], "published_cytology", f"ref_{index}"])
        return path

    def run_builder(self, root: Path, loss: Path, ploidy: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(BUILD),
                "--species-loss-dir",
                str(loss),
                "--ploidy-ledger",
                str(ploidy),
                "--output-dir",
                str(root / "pgls_input"),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_complete_input_and_reports_reconcile(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            loss = self.make_aggregation(root)
            ploidy = self.write_ploidy(root)
            completed = self.run_builder(root, loss, ploidy)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output = root / "pgls_input"
            rows = read_tsv(output / "pgls_input.tsv")
            self.assertEqual([row["biological_species"] for row in rows], SPECIES)
            self.assertEqual(
                [int(row["lineage_specific_nonshared_positive_loss_count"]) for row in rows],
                [0, 1, 2, 3, 4, 5],
            )
            self.assertEqual([int(row["callable_denominator"]) for row in rows], [5] * 6)
            ploidy_report = json.loads((output / "ploidy_ledger_pass.json").read_text())
            self.assertEqual(ploidy_report["status"], "PASS")
            report = json.loads((output / "pgls_input_pass.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["input_data"]["sha256"], digest(output / "pgls_input.tsv"))
            self.assertEqual(
                report["upstream_bindings"]["species_loss_manifest"]["sha256"],
                digest(loss / "species_loss_summary.json"),
            )
            checksums = read_tsv(output / "checksums.sha256.tsv")
            self.assertEqual(len(checksums), 3)

    def test_any_unit_rule_and_ploidy_mismatch_fail_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            loss = self.make_aggregation(root, any_unit=True)
            completed = self.run_builder(root, loss, self.write_ploidy(root))
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("requires all_units_positive", completed.stderr)
            self.assertFalse((root / "pgls_input").exists())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            loss = self.make_aggregation(root)
            wrong_species = SPECIES[:-1] + ["Actinidia wrong"]
            completed = self.run_builder(root, loss, self.write_ploidy(root, wrong_species))
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("ploidy ledger order/set differs", completed.stderr)
            self.assertFalse((root / "pgls_input").exists())


if __name__ == "__main__":
    unittest.main()
