"""Tests for legacy expected-deletion-locus localization."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "prepare_legacy_deleted_locus_spatial_inputs.py"


class PrepareLegacyDeletedSpatialInputsTests(unittest.TestCase):
    def test_localizes_deleted_gene_between_bilateral_anchors(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "reference.coords").write_text(
                "g1\tRef1\t1\t10\t+\ng2\tRef1\t20\t30\t+\ng3\tRef1\t40\t50\t+\n",
                encoding="utf-8",
            )
            pairs = root / "pairs.tsv"
            pairs.write_text(
                "g1\tRef1\t1\t10\tt1\tTargetA\t100\t110\n"
                "g3\tRef1\t40\t50\tt3\tTargetA\t300\t310\n",
                encoding="utf-8",
            )
            (root / "deleted.txt").write_text("g2\n", encoding="utf-8")
            (root / "genome.fa").write_text(">TargetA\n" + "A" * 500 + "\n", encoding="utf-8")
            (root / "genes.gff3").write_text(
                "##gff-version 3\nTargetA\ttest\tgene\t50\t60\t.\t+\t.\tID=x\n",
                encoding="utf-8",
            )
            audit = {
                "schema_version": 2, "sample": "Sample",
                "inputs": {
                    "pairs": str(pairs), "reference_column_1_based": 1,
                    "query_column_1_based": 5,
                },
                "metrics": {"unique_reference_anchors": 2},
            }
            (root / "audit.json").write_text(json.dumps(audit), encoding="utf-8")
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "assembly_unit_id\tlegacy_sample\tbiological_species\thaplotype_or_subgenome\tassembly_scope\tsynorth_audit\tsynorth_pairs\tdeleted_genes\tgenome\tgff\n"
                "u\tSample\tActinidia test\tA\tchromosome\taudit.json\tpairs.tsv\tdeleted.txt\tgenome.fa\tgenes.gff3\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(manifest),
                    "--data-root", str(root), "--reference-coords", "reference.coords",
                    "--output-dir", str(root / "out"), "--padding-bp", "0",
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (root / "out" / "expected_deleted_locus_coordinates.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["reference_gene_id"], "g2")
            self.assertEqual(rows[0]["chromosome"], "TargetA")
            self.assertEqual(rows[0]["expected_locus_start_1based"], "100")
            self.assertEqual(rows[0]["expected_locus_end_1based"], "310")
            report = json.loads((root / "out" / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["spatially_localized_positive_deleted_count"], 1)


if __name__ == "__main__":
    unittest.main()
