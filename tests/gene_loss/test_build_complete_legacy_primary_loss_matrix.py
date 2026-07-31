"""Tests for conservative conversion of exact-bound legacy loss evidence."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "build_complete_legacy_primary_loss_matrix.py"


class CompleteLegacyPrimaryMatrixTests(unittest.TestCase):
    def test_primary_and_historical_semantics_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as legacy:
            root = Path(temporary)
            legacy_root = Path(legacy)
            reference_target = legacy_root / "reference.fa"
            reference_target.write_text(">g1\nATG\n>g2\nATG\n>g3\nATG\n>g4\nATG\n", encoding="utf-8")
            (root / "reference.fa").symlink_to(reference_target)
            pairs_target = legacy_root / "pairs.tsv"
            pairs_target.write_text("g1\tx\ttarget\n", encoding="utf-8")
            (root / "pairs.tsv").symlink_to(pairs_target)
            (root / "decayed.txt").write_text("g2\n", encoding="utf-8")
            (root / "deleted.txt").write_text("g3\n", encoding="utf-8")
            audit = {
                "schema_version": 2,
                "sample": "Sample_A",
                "inputs": {
                    "pairs": str(pairs_target),
                    "pairs_has_header": False,
                    "reference_column_1_based": 1,
                },
                "metrics": {
                    "duplicate_pair_rows": 0,
                    "unique_reference_anchors": 1,
                    "reference_anchor_ids_absent_from_fasta_count": 0,
                    "reference_anchor_ids_absent_from_coordinates_count": 0,
                },
            }
            (root / "audit.json").write_text(json.dumps(audit) + "\n", encoding="utf-8")
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "\t".join(
                    [
                        "assembly_unit_id", "legacy_sample", "synorth_audit", "synorth_pairs",
                        "decayed_genes", "deleted_genes",
                    ]
                )
                + "\n"
                + "u1\tSample_A\taudit.json\tpairs.tsv\tdecayed.txt\tdeleted.txt\n",
                encoding="utf-8",
            )
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(manifest),
                    "--data-root", str(root), "--reference-cds", "reference.fa",
                    "--output-dir", str(output),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "complete_unit_loss_matrix.tsv").open(newline="") as handle:
                primary = {row["reference_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(primary["g1"]["classification"], "retained")
            self.assertEqual(primary["g2"]["classification"], "uncertain")
            self.assertEqual(primary["g2"]["callable"], "true")
            self.assertEqual(primary["g3"]["classification"], "deleted")
            self.assertEqual(primary["g4"]["classification"], "uncertain")
            self.assertEqual(primary["g4"]["callable"], "false")
            with (output / "historical_reproduction_loss_matrix.tsv").open(newline="") as handle:
                historical = {
                    row["reference_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(historical["g2"]["classification"], "pseudogenized")
            self.assertEqual(historical["g2"]["positive_loss"], "true")
            self.assertEqual(historical["g3"]["classification"], "deleted")
            self.assertEqual(json.loads((output / "run_manifest.json").read_text())["status"], "PASS")

    def test_overlap_with_synorth_anchor_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "reference.fa").write_text(">g1\nATG\n", encoding="utf-8")
            (root / "pairs.tsv").write_text("g1\tx\n", encoding="utf-8")
            (root / "decayed.txt").write_text("g1\n", encoding="utf-8")
            (root / "deleted.txt").write_text("", encoding="utf-8")
            audit = {
                "schema_version": 2, "sample": "S",
                "inputs": {"pairs": str(root / "pairs.tsv"), "pairs_has_header": False, "reference_column_1_based": 1},
                "metrics": {
                    "duplicate_pair_rows": 0, "unique_reference_anchors": 1,
                    "reference_anchor_ids_absent_from_fasta_count": 0,
                    "reference_anchor_ids_absent_from_coordinates_count": 0,
                },
            }
            (root / "audit.json").write_text(json.dumps(audit))
            (root / "manifest.tsv").write_text(
                "assembly_unit_id\tlegacy_sample\tsynorth_audit\tsynorth_pairs\tdecayed_genes\tdeleted_genes\n"
                "u\tS\taudit.json\tpairs.tsv\tdecayed.txt\tdeleted.txt\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(root / "manifest.tsv"),
                    "--data-root", str(root), "--reference-cds", "reference.fa",
                    "--output-dir", str(root / "out"),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("overlap SynOrths anchors", completed.stderr)


if __name__ == "__main__":
    unittest.main()
