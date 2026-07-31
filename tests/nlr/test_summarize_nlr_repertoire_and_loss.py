"""Tests for exact-cohort NLR repertoire and positive-loss summaries."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY_ROOT / "scripts" / "nlr" / "summarize_nlr_repertoire_and_loss.py"
DIGEST_A = "a" * 64
DIGEST_B = "b" * 64


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SummarizeNlrRepertoireAndLossTest(unittest.TestCase):
    def write_fixture(self, root: Path) -> dict[str, Path]:
        metadata = root / "assembly_units.tsv"
        repertoire = root / "complete_nlr_counts.tsv"
        calls = root / "positive_reference_nlr_calls.tsv"
        denominators = root / "callable_reference_nlr_denominators.tsv"
        metadata.write_text(
            "assembly_unit_id\tbiological_species\thaplotype_or_subgenome\tassembly_scope\tinclude\tanalysis_cohort\n"
            "deliciosa_a\tActinidia deliciosa\tA\tchromosome_partition\ttrue\tprimary20\n"
            "deliciosa_b\tActinidia deliciosa\tB\tchromosome_partition\ttrue\tprimary20\n"
            "rufa_chr\tActinidia rufa\tunphased\tchromosome_only_subset\ttrue\trufa_sensitivity21\n",
            encoding="utf-8",
        )
        repertoire.write_text(
            "assembly_unit_id\tassembly_scope\ttotal_nlr_count\trepertoire_source_basename\trepertoire_source_sha256\n"
            f"deliciosa_a\tchromosome_partition\t100\tdeliciosa_a_nlr.tsv\t{DIGEST_A}\n"
            f"deliciosa_b\tchromosome_partition\t90\tdeliciosa_b_nlr.tsv\t{DIGEST_B}\n",
            encoding="utf-8",
        )
        calls.write_text(
            "assembly_unit_id\tassembly_scope\treference_nlr_id\treference_nlr_universe_id\n"
            "deliciosa_a\tchromosome_partition\tREF_NLR_1\tcs_complete_nlr_v1\n"
            "deliciosa_a\tchromosome_partition\tREF_NLR_2\tcs_complete_nlr_v1\n",
            encoding="utf-8",
        )
        denominators.write_text(
            "assembly_unit_id\tassembly_scope\tcallable_reference_nlr_denominator\t"
            "reference_nlr_universe_id\tdenominator_source_basename\tdenominator_source_sha256\n"
            f"deliciosa_a\tchromosome_partition\t10\tcs_complete_nlr_v1\tcallable_a.tsv\t{DIGEST_A}\n"
            f"deliciosa_b\tchromosome_partition\t0\tcs_complete_nlr_v1\tcallable_b.tsv\t{DIGEST_B}\n",
            encoding="utf-8",
        )
        return {
            "metadata": metadata,
            "repertoire": repertoire,
            "calls": calls,
            "denominators": denominators,
        }

    def run_script(self, root: Path, fixture: dict[str, Path], *, output_name: str = "summary"):
        output = root / output_name
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--metadata",
                str(fixture["metadata"]),
                "--repertoire-counts",
                str(fixture["repertoire"]),
                "--positive-loss-calls",
                str(fixture["calls"]),
                "--callable-denominators",
                str(fixture["denominators"]),
                "--analysis-cohort",
                "primary20",
                "--cohort-role",
                "primary",
                "--output-dir",
                str(output),
            ],
            cwd=REPOSITORY_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        return completed, output

    def test_count_denominators_zero_calls_and_species_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self.write_fixture(root)
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {
                    "nlr_unit_summary.tsv",
                    "nlr_species_aggregate.tsv",
                    "validation.json",
                    "input_checksums.tsv",
                },
            )

            units = read_tsv(output / "nlr_unit_summary.tsv")
            self.assertEqual([row["assembly_unit_id"] for row in units], ["deliciosa_a", "deliciosa_b"])
            self.assertEqual(units[0]["total_nlr_count"], "100")
            self.assertEqual(units[0]["positive_reference_nlr_loss_count"], "2")
            self.assertEqual(float(units[0]["positive_reference_nlr_loss_percentage"]), 20.0)
            self.assertEqual(units[1]["positive_reference_nlr_loss_count"], "0")
            self.assertEqual(units[1]["callable_reference_nlr_denominator"], "0")
            self.assertEqual(units[1]["positive_reference_nlr_loss_percentage"], "")
            self.assertEqual(units[1]["percentage_status"], "undefined_zero_denominator")
            self.assertNotIn("terminal", units[0])
            self.assertNotIn("n_lost_genes", units[0])

            species = read_tsv(output / "nlr_species_aggregate.tsv")
            self.assertEqual(len(species), 1)
            self.assertEqual(species[0]["assembly_unit_count"], "2")
            self.assertEqual(species[0]["total_nlr_count_sum_across_units"], "190")
            self.assertEqual(species[0]["positive_reference_nlr_loss_count_sum_across_units"], "2")
            self.assertEqual(float(species[0]["positive_reference_nlr_loss_percentage_across_unit_comparisons"]), 20.0)

            validation = json.loads((output / "validation.json").read_text(encoding="utf-8"))
            self.assertEqual(validation["denominator_input_mode"], "count")
            self.assertEqual(validation["undefined_percentage_unit_count"], 1)
            checksums_text = (output / "input_checksums.tsv").read_text(encoding="utf-8")
            self.assertNotIn(str(root), checksums_text)

    def test_positive_count_cannot_exceed_count_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self.write_fixture(root)
            fixture["denominators"].write_text(
                "assembly_unit_id\tassembly_scope\tcallable_reference_nlr_denominator\t"
                "reference_nlr_universe_id\tdenominator_source_basename\tdenominator_source_sha256\n"
                f"deliciosa_a\tchromosome_partition\t1\tcs_complete_nlr_v1\tcallable_a.tsv\t{DIGEST_A}\n"
                f"deliciosa_b\tchromosome_partition\t0\tcs_complete_nlr_v1\tcallable_b.tsv\t{DIGEST_B}\n",
                encoding="utf-8",
            )
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("exceeds callable denominator", completed.stderr)
            self.assertFalse(output.exists())

    def test_catalog_denominator_rejects_positive_call_outside_callable_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self.write_fixture(root)
            fixture["denominators"].write_text(
                "assembly_unit_id\tassembly_scope\treference_nlr_id\treference_nlr_universe_id\t"
                "denominator_source_basename\tdenominator_source_sha256\n"
                f"deliciosa_a\tchromosome_partition\tREF_NLR_1\tcs_complete_nlr_v1\tcallable_a.tsv\t{DIGEST_A}\n"
                f"deliciosa_a\tchromosome_partition\tREF_NLR_3\tcs_complete_nlr_v1\tcallable_a.tsv\t{DIGEST_A}\n"
                f"deliciosa_b\tchromosome_partition\tREF_NLR_1\tcs_complete_nlr_v1\tcallable_b.tsv\t{DIGEST_B}\n",
                encoding="utf-8",
            )
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("absent from the callable denominator", completed.stderr)
            self.assertFalse(output.exists())

    def test_exact_cohort_and_unique_metadata_rows_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self.write_fixture(root)
            with fixture["repertoire"].open("a", encoding="utf-8") as handle:
                handle.write(
                    f"rufa_chr\tchromosome_only_subset\t50\trufa.tsv\t{DIGEST_A}\n"
                )
            completed, output = self.run_script(root, fixture, output_name="extra")
            self.assertEqual(completed.returncode, 2)
            self.assertIn("does not exactly match the included cohort", completed.stderr)
            self.assertFalse(output.exists())

            fixture = self.write_fixture(root)
            with fixture["metadata"].open("a", encoding="utf-8") as handle:
                handle.write(
                    "deliciosa_a\tActinidia deliciosa\tA\tchromosome_partition\ttrue\tprimary20\n"
                )
            completed, output = self.run_script(root, fixture, output_name="duplicate")
            self.assertEqual(completed.returncode, 2)
            self.assertIn("duplicate assembly_unit_id", completed.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
