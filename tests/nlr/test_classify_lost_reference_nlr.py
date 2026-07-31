from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "nlr" / "classify_lost_reference_nlr.py"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class LostReferenceNlrClassTests(unittest.TestCase):
    def test_structural_class_counts_and_denominators_close(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.tsv"
            universe = root / "universe.tsv"
            calls = root / "calls.txt"
            matrix = root / "matrix.tsv"
            nlr_output_root = root / "nlr"
            unit_summary = root / "unit_summary.tsv"
            output = root / "out"
            write_tsv(
                metadata,
                [
                    "assembly_unit_id",
                    "biological_species",
                    "haplotype_or_subgenome",
                    "include",
                ],
                [
                    {
                        "assembly_unit_id": "u1",
                        "biological_species": "Actinidia alpha",
                        "haplotype_or_subgenome": "A",
                        "include": "true",
                    },
                    {
                        "assembly_unit_id": "u2",
                        "biological_species": "Actinidia beta",
                        "haplotype_or_subgenome": "",
                        "include": "true",
                    },
                ],
            )
            write_tsv(
                universe,
                [
                    "reference_nlr_id",
                    "included_in_article_nonshared_analysis",
                    "exclusion_reason",
                ],
                [
                    {
                        "reference_nlr_id": "g1",
                        "included_in_article_nonshared_analysis": "true",
                        "exclusion_reason": "",
                    },
                    {
                        "reference_nlr_id": "g2",
                        "included_in_article_nonshared_analysis": "true",
                        "exclusion_reason": "",
                    },
                    {
                        "reference_nlr_id": "g3",
                        "included_in_article_nonshared_analysis": "false",
                        "exclusion_reason": (
                            "positive_in_all_23_units_under_article_method"
                        ),
                    },
                ],
            )
            calls.write_text(
                "g1\tg1_nlr1\tCC-NBARC\t1\t10\t+\tmotif_1\n"
                "g2\tg2_nlr1\tTIR-NBARC-LRR\t1\t10\t+\tmotif_1\n"
                "g3\tg3_nlr1\tNBARC\t1\t10\t+\tmotif_1\n",
                encoding="utf-8",
            )
            for unit in ("u1", "u2"):
                (nlr_output_root / unit).mkdir(parents=True)
            (nlr_output_root / "u1" / "nlr_calls.txt").write_text(
                "chr1\tchr1_nlr1\tCC-NBARC\t1\t10\t+\tmotif_1\n"
                "chr1\tchr1_nlr2\tTIR-NBARC-LRR\t20\t30\t+\tmotif_1\n",
                encoding="utf-8",
            )
            (nlr_output_root / "u2" / "nlr_calls.txt").write_text(
                "chr2\tchr2_nlr1\tNBARC\t1\t10\t+\tmotif_1\n",
                encoding="utf-8",
            )
            write_tsv(
                unit_summary,
                ["assembly_unit_id", "total_nlr_count"],
                [
                    {"assembly_unit_id": "u1", "total_nlr_count": 2},
                    {"assembly_unit_id": "u2", "total_nlr_count": 1},
                ],
            )
            matrix_rows = [
                {
                    "reference_gene_id": "g1",
                    "assembly_unit_id": "u1",
                    "manuscript_classification": "decayed",
                    "refined_decayed_cause": "local_sequence_no_explicit_coding_disruption",
                    "refined_cause_evidence_level": "sequence_detected_only",
                },
                {
                    "reference_gene_id": "g2",
                    "assembly_unit_id": "u1",
                    "manuscript_classification": "retained",
                    "refined_decayed_cause": "not_applicable_retained",
                    "refined_cause_evidence_level": "synorth_anchor",
                },
                {
                    "reference_gene_id": "g1",
                    "assembly_unit_id": "u2",
                    "manuscript_classification": "deleted",
                    "refined_decayed_cause": "no_qualifying_genomewide_tblastx_hit",
                    "refined_cause_evidence_level": "sequence_absent",
                },
                {
                    "reference_gene_id": "g2",
                    "assembly_unit_id": "u2",
                    "manuscript_classification": "not_called_loss",
                    "refined_decayed_cause": "not_called",
                    "refined_cause_evidence_level": "not_called",
                },
            ]
            write_tsv(matrix, list(matrix_rows[0]), matrix_rows)
            command = [
                sys.executable,
                str(SCRIPT),
                "--loss-matrix",
                str(matrix),
                "--reference-nlr-universe",
                str(universe),
                "--reference-nlr-calls",
                str(calls),
                "--unit-metadata",
                str(metadata),
                "--nlr-output-root",
                str(nlr_output_root),
                "--nlr-unit-summary",
                str(unit_summary),
                "--output-dir",
                str(output),
                "--expected-units",
                "2",
                "--expected-reference-nlrs",
                "3",
                "--expected-nonshared-reference-nlrs",
                "2",
                "--expected-positive-calls",
                "2",
                "--expected-resolved-denominator",
                "3",
                "--expected-repertoire-nlrs",
                "3",
            ]
            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "lost_nlr_structural_class_summary.tsv").open(
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 18)
            indexed = {
                (row["assembly_unit_id"], row["reference_nlr_class"]): row
                for row in rows
            }
            self.assertEqual(indexed[("u1", "CC-NBARC")]["positive_loss_count"], "1")
            self.assertEqual(
                indexed[("u2", "TIR-NBARC-LRR")][
                    "resolved_unit_gene_denominator"
                ],
                "0",
            )
            with (output / "lost_nlr_structural_class_calls.tsv").open(
                encoding="utf-8",
                newline="",
            ) as handle:
                positives = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(positives), 2)
            self.assertEqual(
                {row["reference_nlr_class"] for row in positives},
                {"CC-NBARC"},
            )
            with (
                output / "shared_nlr_structural_class_summary.tsv"
            ).open(encoding="utf-8", newline="") as handle:
                shared = list(csv.DictReader(handle, delimiter="\t"))
            shared_index = {
                row["reference_nlr_class"]: row[
                    "shared_reference_nlr_gene_count"
                ]
                for row in shared
            }
            self.assertEqual(shared_index["NBARC"], "1")
            with (
                output / "nlr_repertoire_structural_class_summary.tsv"
            ).open(encoding="utf-8", newline="") as handle:
                repertoire = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(repertoire), 18)
            self.assertEqual(
                sum(int(row["nlr_gene_count"]) for row in repertoire),
                3,
            )
            with (
                output / "nlr_structural_class_loss_rates.tsv"
            ).open(encoding="utf-8", newline="") as handle:
                rates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rates), 18)
            rate_index = {
                (row["assembly_unit_id"], row["reference_nlr_class"]): row
                for row in rates
            }
            self.assertEqual(
                rate_index[("u1", "CC-NBARC")]["all_loss_percentage"],
                "100.000000",
            )
            self.assertEqual(
                rate_index[("u1", "NBARC")]["all_loss_percentage"],
                "100.000000",
            )
            self.assertEqual(
                rate_index[("u1", "NBARC")]["nonshared_loss_percentage"],
                "",
            )


if __name__ == "__main__":
    unittest.main()
