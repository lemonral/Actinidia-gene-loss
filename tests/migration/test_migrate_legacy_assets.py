"""Tests for size-aware legacy asset migration."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).parents[2] / "scripts" / "migration" / "migrate_legacy_assets.py"


class LegacyMigrationTest(unittest.TestCase):
    def test_small_assets_are_copied_and_large_assets_are_linked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            sources = root / "sources"
            sources.mkdir()
            genome = sources / "sample.fa"
            gff = sources / "sample.gff3"
            protein = sources / "sample.faa"
            genome.write_bytes(b"G" * 100)
            gff.write_bytes(b"A" * 20)
            protein.write_bytes(b"P" * 10)

            mapping = root / "mapping.tsv"
            mapping.write_text(
                "legacy_sample\tassembly_unit_id\tbiological_species\t"
                "haplotype_or_subgenome\taccession\tsource_url\n"
                "Old_A\tact_species_a\tActinidia species\tA\tACC1\thttps://example.org/ACC1\n",
                encoding="utf-8",
            )
            legacy = root / "legacy.tsv"
            legacy.write_text(
                "sample\tgenome\tgff\tprotein\n"
                f"Old_A\t{genome}\t{gff}\t{protein}\n",
                encoding="utf-8",
            )
            data_root = root / "data"
            report = root / "report.json"
            resolved_manifest = root / "legacy_resolved.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--mapping",
                    str(mapping),
                    "--legacy-manifest",
                    str(legacy),
                    "--data-root",
                    str(data_root),
                    "--report",
                    str(report),
                    "--resolved-manifest",
                    str(resolved_manifest),
                    "--copy-max-mib",
                    "0.00005",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            migrated = data_root / "legacy_linked" / "act_species_a"
            self.assertTrue((migrated / "genome.fa").is_symlink())
            self.assertFalse((migrated / "gff.gff3").is_symlink())
            self.assertFalse((migrated / "protein.faa").is_symlink())
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_count"], 1)
            self.assertEqual(payload["asset_count"], 3)
            self.assertEqual({row["mode"] for row in payload["assets"]}, {"copy", "symlink"})
            self.assertTrue(all(len(row["sha256"]) == 64 for row in payload["assets"]))
            header, row = resolved_manifest.read_text(encoding="utf-8").splitlines()
            self.assertIn("target_haplotype", header.split("\t"))
            self.assertIn("include_downstream", header.split("\t"))
            self.assertIn("act_species_a", row.split("\t"))

    def test_sample_set_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            mapping.write_text(
                "legacy_sample\tassembly_unit_id\tbiological_species\t"
                "haplotype_or_subgenome\taccession\tsource_url\n"
                "Old_A\tact_species_a\tActinidia species\tA\tACC1\thttps://example.org/ACC1\n",
                encoding="utf-8",
            )
            legacy = root / "legacy.tsv"
            legacy.write_text("sample\tgenome\tgff\tprotein\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--mapping",
                    str(mapping),
                    "--legacy-manifest",
                    str(legacy),
                    "--data-root",
                    str(root / "data"),
                    "--report",
                    str(root / "report.json"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("sample sets differ", completed.stderr)


if __name__ == "__main__":
    unittest.main()
