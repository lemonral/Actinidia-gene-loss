from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "scripts"
    / "function"
    / "prepare_unit_article_method_foregrounds.py"
)
SPEC = importlib.util.spec_from_file_location("unit_foregrounds", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class UnitArticleMethodForegroundTests(unittest.TestCase):
    def test_units_are_independent_and_not_called_is_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proteins = root / "reference.faa"
            proteins.write_text(
                ">g1\nM\n>g2\nM\n>g3\nM\n",
                encoding="utf-8",
            )
            metadata = root / "units.tsv"
            metadata.write_text(
                "assembly_unit_id\tbiological_species\t"
                "haplotype_or_subgenome\tinclude\n"
                "u1\tActinidia alpha\tA\ttrue\n"
                "u2\tActinidia beta\tB\ttrue\n",
                encoding="utf-8",
            )
            matrix = root / "matrix.tsv.gz"
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "reference_gene_id",
                        "assembly_unit_id",
                        "manuscript_classification",
                        "manuscript_positive_loss",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for unit, states in (
                    ("u1", ("retained", "decayed", "not_called_loss")),
                    ("u2", ("deleted", "retained", "decayed")),
                ):
                    for gene, state in zip(("g1", "g2", "g3"), states):
                        writer.writerow(
                            {
                                "reference_gene_id": gene,
                                "assembly_unit_id": unit,
                                "manuscript_classification": state,
                                "manuscript_positive_loss": str(
                                    state in {"decayed", "deleted"}
                                ).lower(),
                            }
                        )
            output = root / "out"
            MODULE.run(
                argparse.Namespace(
                    unit_matrix=matrix,
                    unit_metadata=metadata,
                    reference_protein=proteins,
                    expected_units=2,
                    expected_reference_genes=3,
                    expected_matrix_sha256="",
                    output_dir=output,
                )
            )
            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["status"],
                "PASS_UNIT_ARTICLE_METHOD_FOREGROUNDS",
            )
            self.assertEqual(
                manifest["counts"]["foreground_memberships"],
                3,
            )
            self.assertEqual(
                manifest["counts"]["resolved_background_memberships"],
                5,
            )
            with gzip.open(
                output / "foreground_gene_ids.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                foregrounds = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {
                    (row["foreground_id"], row["reference_gene_id"])
                    for row in foregrounds
                },
                {("unit__u1", "g2"), ("unit__u2", "g1"), ("unit__u2", "g3")},
            )


if __name__ == "__main__":
    unittest.main()
