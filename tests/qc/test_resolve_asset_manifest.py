"""Fail-closed tests for the canonical analysis-unit manifest generator."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "resolve_asset_manifest.py"
SPEC = importlib.util.spec_from_file_location("resolve_asset_manifest", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

ASSEMBLY_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "ploidy",
    "accession",
    "version",
    "source_bundle_id",
    "assembly_scope",
    "partition_rule",
    "genome_url",
    "annotation_url",
    "protein_url",
    "expected_genome_sha256",
    "expected_annotation_sha256",
    "expected_protein_sha256",
    "qc_status",
    "include_qc",
    "include_gene_loss",
    "include_species_tree",
    "exclusion_reason",
    "notes",
)
DOWNLOAD_COLUMNS = (
    "asset_id",
    "assembly_unit_id",
    "asset_type",
    "url",
    "relative_path",
    "expected_bytes",
    "md5",
    "sha256",
    "download",
    "source_note",
)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ResolutionFixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.data_root = root / "data"
        self.data_root.mkdir()
        self.assemblies = root / "assemblies.tsv"
        self.downloads = root / "downloads.tsv"
        self.report = root / "download_report.json"
        self.output = root / "resolved.tsv"

        self.payloads = {
            "genome": b">chr1\nACGTNN\n",
            "gff": b"##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
            "protein": b">p1\nMAAA\n",
        }
        self.urls = {
            role: f"https://example.org/unit.{role}" for role in self.payloads
        }
        self.assembly_row = {
            "assembly_unit_id": "act_test_hap1",
            "biological_species": "Actinidia testensis",
            "individual_id": "individual_1",
            "haplotype_or_subgenome": "HAP1",
            "ploidy": "2x",
            "accession": "GCA_000000001.1",
            "version": "1",
            "source_bundle_id": "test_bundle",
            "assembly_scope": "complete_test_scope",
            "partition_rule": "none",
            "genome_url": self.urls["genome"],
            "annotation_url": self.urls["gff"],
            "protein_url": self.urls["protein"],
            "expected_genome_sha256": hashlib.sha256(self.payloads["genome"]).hexdigest(),
            "expected_annotation_sha256": "",
            "expected_protein_sha256": "",
            "qc_status": "ready",
            "include_qc": "true",
            "include_gene_loss": "true",
            "include_species_tree": "false",
            "exclusion_reason": "",
            "notes": "test row",
        }
        self.download_rows: list[dict[str, str]] = []
        self.report_rows: list[dict[str, object]] = []
        for role, payload in self.payloads.items():
            relative = f"downloads/test/unit.{role}"
            path = self.data_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            md5 = hashlib.md5(payload).hexdigest() if role == "genome" else ""
            sha256 = hashlib.sha256(payload).hexdigest()
            self.download_rows.append(
                {
                    "asset_id": f"test_{role}",
                    "assembly_unit_id": "act_test_hap1",
                    "asset_type": role,
                    "url": self.urls[role],
                    "relative_path": relative,
                    "expected_bytes": str(len(payload)),
                    "md5": md5,
                    "sha256": "",
                    "download": "true",
                    "source_note": "test asset",
                }
            )
            self.report_rows.append(
                {
                    "asset_id": f"test_{role}",
                    "assembly_unit_id": "act_test_hap1",
                    "asset_type": role,
                    "relative_path": relative,
                    "bytes": len(payload),
                    "md5": md5,
                    "sha256": sha256,
                    "publisher_checksum_declared": bool(md5),
                    "verified": bool(md5),
                    "route": "direct",
                    "status": "verified",
                    "problems": [],
                }
            )

    def write(self) -> None:
        write_tsv(self.assemblies, ASSEMBLY_COLUMNS, [self.assembly_row])
        write_tsv(self.downloads, DOWNLOAD_COLUMNS, self.download_rows)
        self.report.write_text(json.dumps(self.report_rows), encoding="utf-8")

    def resolve(self) -> list[dict[str, str]]:
        return MODULE.resolve_manifests(
            assemblies_path=self.assemblies,
            downloads_path=self.downloads,
            data_root=self.data_root,
            report_path=self.report,
            output_path=self.output,
        )


class ResolveAssetManifestTest(unittest.TestCase):
    def test_valid_assets_emit_canonical_and_compatibility_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.write()
            rows = fixture.resolve()

            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["assembly_unit_id"], "act_test_hap1")
            self.assertEqual(row["sample"], "act_test_hap1")
            self.assertEqual(row["target_haplotype"], "act_test_hap1")
            self.assertEqual(row["species"], "Actinidia testensis")
            self.assertEqual(row["include_downstream"], "true")
            self.assertEqual(row["current_or_alternative"], "current")
            self.assertEqual(row["genome_integrity_level"], "publisher_md5")
            self.assertEqual(row["gff_integrity_level"], "expected_size_plus_local_sha256")
            self.assertEqual(Path(row["protein"]), (fixture.data_root / "downloads/test/unit.protein").resolve())

            with fixture.output.open(encoding="utf-8", newline="") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(written, rows)
            self.assertEqual(tuple(written[0]), MODULE.OUTPUT_COLUMNS)

    def test_duplicate_enabled_role_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            duplicate = dict(fixture.download_rows[0])
            duplicate["asset_id"] = "duplicate_genome"
            duplicate["relative_path"] = "downloads/test/duplicate.genome"
            (fixture.data_root / duplicate["relative_path"]).write_bytes(fixture.payloads["genome"])
            fixture.download_rows.append(duplicate)
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "duplicate enabled genome role"):
                fixture.resolve()

    def test_missing_required_role_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.download_rows = [row for row in fixture.download_rows if row["asset_type"] != "protein"]
            fixture.report_rows = [row for row in fixture.report_rows if row["asset_type"] != "protein"]
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "missing enabled required asset roles: protein"):
                fixture.resolve()

    def test_local_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.write()
            (fixture.data_root / "downloads/test/unit.gff").write_bytes(b"changed")
            with self.assertRaisesRegex(MODULE.ResolutionError, "recorded bytes|local SHA-256"):
                fixture.resolve()

    def test_publisher_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.download_rows[0]["md5"] = "0" * 32
            fixture.report_rows[0]["md5"] = "0" * 32
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "publisher MD5"):
                fixture.resolve()

    def test_incomplete_report_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.report_rows.pop()
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "does not exactly cover enabled assets"):
                fixture.resolve()

    def test_parent_path_escape_is_rejected_even_for_disabled_asset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.download_rows.append(
                {
                    "asset_id": "disabled_escape",
                    "assembly_unit_id": "act_test_hap1",
                    "asset_type": "cds",
                    "url": "https://example.org/escape",
                    "relative_path": "../escape.fa",
                    "expected_bytes": "1",
                    "md5": "",
                    "sha256": "",
                    "download": "false",
                    "source_note": "must still be safe",
                }
            )
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "escapes --data-root"):
                fixture.resolve()

    def test_failed_downloader_record_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = ResolutionFixture(Path(temporary_directory))
            fixture.report_rows.append({"status": "failed", "error": "transfer stopped"})
            fixture.write()
            with self.assertRaisesRegex(MODULE.ResolutionError, "downloader recorded failure"):
                fixture.resolve()


if __name__ == "__main__":
    unittest.main()
