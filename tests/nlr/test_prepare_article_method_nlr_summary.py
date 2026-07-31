"""Regression tests for the article-method NLR adapter."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "nlr" / "prepare_article_method_nlr_summary.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


class ArticleMethodNlrSummaryTest(unittest.TestCase):
    def test_builds_nonshared_resolved_unit_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.tsv"
            write_tsv(
                metadata,
                [
                    {
                        "assembly_unit_id": "u1",
                        "biological_species": "Actinidia alpha",
                        "haplotype_or_subgenome": "A",
                        "assembly_scope": "scope",
                        "include": "true",
                    },
                    {
                        "assembly_unit_id": "u2",
                        "biological_species": "Actinidia beta",
                        "haplotype_or_subgenome": "unphased",
                        "assembly_scope": "scope",
                        "include": "true",
                    },
                ],
            )
            universe = root / "universe.tsv"
            write_tsv(
                universe,
                [
                    {"reference_nlr_id": "g1"},
                    {"reference_nlr_id": "g2"},
                ],
            )
            shared = root / "shared.tsv"
            write_tsv(shared, [{"reference_gene_id": "g1"}])
            repertoire = root / "repertoire.tsv"
            write_tsv(
                repertoire,
                [
                    {"assembly_unit_id": "u1", "total_nlr_count": "10"},
                    {"assembly_unit_id": "u2", "total_nlr_count": "11"},
                ],
            )
            matrix = root / "matrix.tsv.gz"
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "reference_gene_id",
                        "assembly_unit_id",
                        "manuscript_classification",
                        "refined_decayed_cause",
                        "refined_cause_evidence_level",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "reference_gene_id": gene,
                            "assembly_unit_id": unit,
                            "manuscript_classification": classification,
                            "refined_decayed_cause": (
                                "not_applicable_retained"
                                if classification == "retained"
                                else (
                                    "frameshift_supported"
                                    if classification == "decayed"
                                    else "no_qualifying_genomewide_tblastx_hit"
                                )
                            ),
                            "refined_cause_evidence_level": (
                                "exact_synorth"
                                if classification == "retained"
                                else (
                                    "explicit_coding_disruption"
                                    if classification == "decayed"
                                    else "manuscript_threshold"
                                )
                            ),
                        }
                        for unit, classes in (
                            ("u1", ("decayed", "decayed", "retained")),
                            ("u2", ("deleted", "retained", "retained")),
                        )
                        for gene, classification in zip(
                            ("g1", "g2", "g3"),
                            classes,
                        )
                    ]
                )
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--article-matrix",
                    str(matrix),
                    "--shared-genes",
                    str(shared),
                    "--reference-nlr-universe",
                    str(universe),
                    "--repertoire-counts",
                    str(repertoire),
                    "--unit-metadata",
                    str(metadata),
                    "--output-dir",
                    str(output),
                    "--expected-units",
                    "2",
                    "--expected-species",
                    "2",
                    "--expected-reference-genes",
                    "3",
                    "--expected-reference-nlrs",
                    "2",
                    "--expected-nonshared-reference-nlrs",
                    "1",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "PASS_ARTICLE_METHOD_NLR_SUMMARY")
            self.assertEqual(manifest["article_shared_reference_nlrs_excluded"], 1)
            self.assertEqual(manifest["positive_unit_gene_calls"], 1)
            with (output / "article_nlr_unit_summary.tsv").open(
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["article_decayed_reference_nlr_loss_count"], "1")
            self.assertEqual(rows[0]["frameshift_supported_count"], "1")
            self.assertEqual(rows[1]["article_retained_reference_nlr_count"], "1")
            with (output / "nlr_loss_type_summary.tsv").open(
                encoding="utf-8",
                newline="",
            ) as handle:
                type_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(type_rows), 12)
            self.assertEqual(
                sum(int(row["positive_reference_nlr_loss_count"]) for row in type_rows),
                1,
            )


if __name__ == "__main__":
    unittest.main()
