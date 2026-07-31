from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "migration" / "migrate_reference_assets.py"
SPEC = importlib.util.spec_from_file_location("migrate_reference_assets", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ReferenceMigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.legacy = self.root / "legacy"
        self.data = self.root / "data"
        self.legacy.mkdir()
        self.data.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def write_mapping(self, source: Path, checksum: str, destination: str = "legacy_linked/ref/x.fa") -> Path:
        mapping = self.root / "mapping.tsv"
        relative = source.relative_to(self.legacy).as_posix()
        mapping.write_text(
            "reference_id\tbiological_species\trole\tlegacy_relative_path\t"
            "data_relative_path\texpected_sha256\tstatus\tnotes\n"
            f"ref\tSpecies name\tgenome\t{relative}\t{destination}\t{checksum}\tpending\ttest\n",
            encoding="utf-8",
        )
        return mapping

    def args(self, mapping: Path, copy_max_mib: float) -> argparse.Namespace:
        return argparse.Namespace(
            mapping=mapping,
            legacy_root=self.legacy,
            data_root=self.data,
            report=self.root / "report.json",
            resolved_manifest=self.root / "resolved.tsv",
            copy_max_mib=copy_max_mib,
            dry_run=False,
        )

    def test_small_asset_is_copied_and_reconciled(self):
        source = self.legacy / "source.fa"
        source.write_bytes(b">x\nACGT\n")
        checksum = hashlib.sha256(source.read_bytes()).hexdigest()
        payload = MODULE.run(self.args(self.write_mapping(source, checksum), 1.0))
        destination = self.data / "legacy_linked/ref/x.fa"
        self.assertTrue(destination.is_file())
        self.assertFalse(destination.is_symlink())
        self.assertEqual(payload["assets"][0]["mode"], "copy")
        self.assertTrue((self.root / "resolved.tsv").is_file())

    def test_zero_copy_limit_creates_expected_soft_link(self):
        source = self.legacy / "source.fa"
        source.write_bytes(b">x\nACGT\n")
        checksum = hashlib.sha256(source.read_bytes()).hexdigest()
        MODULE.run(self.args(self.write_mapping(source, checksum), 0.0))
        destination = self.data / "legacy_linked/ref/x.fa"
        self.assertTrue(destination.is_symlink())
        self.assertEqual(destination.resolve(), source.resolve())
        # A second run must validate and retain the existing link rather than
        # following it during the data-root containment check.
        MODULE.run(self.args(self.write_mapping(source, checksum), 0.0))
        self.assertTrue(destination.is_symlink())

    def test_checksum_mismatch_fails_before_destination(self):
        source = self.legacy / "source.fa"
        source.write_bytes(b">x\nACGT\n")
        with self.assertRaises(MODULE.ReferenceMigrationError):
            MODULE.run(self.args(self.write_mapping(source, "0" * 64), 1.0))
        self.assertFalse((self.data / "legacy_linked/ref/x.fa").exists())

    def test_parent_path_is_rejected(self):
        source = self.legacy / "source.fa"
        source.write_bytes(b">x\nACGT\n")
        checksum = hashlib.sha256(source.read_bytes()).hexdigest()
        mapping = self.write_mapping(source, checksum, "../escape.fa")
        with self.assertRaises(MODULE.ReferenceMigrationError):
            MODULE.run(self.args(mapping, 1.0))


if __name__ == "__main__":
    unittest.main()
