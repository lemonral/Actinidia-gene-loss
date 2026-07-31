"""Contracts for exact public rooting/calibration-evaluation outgroup assets."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "phylogeny"
MANIFEST_COLUMNS = [
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
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


class PublicOutgroupManifestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.header, cls.rows = read_tsv(CONFIG / "public_outgroup_downloads.tsv")
        _, taxa = read_tsv(CONFIG / "taxa.tsv")
        cls.taxa = {row["taxon_id"]: row for row in taxa}
        _, brackets = read_tsv(CONFIG / "fossil_bracketing_taxa.tsv")
        cls.brackets = {row["taxon_or_lineage"]: row for row in brackets}

    def test_manifest_schema_identity_and_integrity_fields(self) -> None:
        self.assertEqual(self.header, MANIFEST_COLUMNS)
        self.assertEqual(len(self.rows), 10)
        self.assertEqual(
            len({row["asset_id"] for row in self.rows}), len(self.rows)
        )
        self.assertEqual(
            len({row["relative_path"] for row in self.rows}), len(self.rows)
        )
        for row in self.rows:
            self.assertEqual(row["download"], "true")
            self.assertGreater(int(row["expected_bytes"]), 0)
            self.assertEqual(len(row["md5"]), 32)
            self.assertTrue(all(character in "0123456789abcdef" for character in row["md5"]))
            self.assertIn(urlparse(row["url"]).scheme, {"http", "https"})
            self.assertFalse(Path(row["relative_path"]).is_absolute())
            self.assertNotIn("..", Path(row["relative_path"]).parts)

    def test_exact_genome_derived_outgroup_bundles_are_frozen(self) -> None:
        by_unit: dict[str, list[dict[str, str]]] = {}
        for row in self.rows:
            by_unit.setdefault(row["assembly_unit_id"], []).append(row)

        leea = by_unit["leea_coccinea_1464"]
        self.assertEqual({row["asset_type"] for row in leea}, {"genome", "gff"})
        self.assertTrue(all("zenodo.org/api/records/13362874/" in row["url"] for row in leea))
        self.assertEqual({row["assembly_unit_id"] for row in leea}, {"leea_coccinea_1464"})
        self.assertIn("one declared representative haplotype", leea[0]["source_note"])

        catharanthus = by_unit["catharanthus_roseus_asm2450571"]
        self.assertEqual(
            {row["asset_type"] for row in catharanthus},
            {
                "genome",
                "gff",
                "protein",
                "cds",
                "assembly_report",
                "assembly_stats",
            },
        )
        self.assertTrue(
            all("GCA_024505715.1_ASM2450571v1" in row["url"] for row in catharanthus)
        )

    def test_saurauia_is_tree_only_and_excluded_from_count_matrices(self) -> None:
        saurauia_assets = [
            row
            for row in self.rows
            if row["assembly_unit_id"] == "saurauia_tristyla_tree_only"
        ]
        self.assertEqual(
            {row["asset_type"] for row in saurauia_assets},
            {"published_tree_archive", "metadata"},
        )
        self.assertTrue(
            all("downloads/phylogeny_tree_only/" in row["relative_path"] for row in saurauia_assets)
        )

        taxon = self.taxa["saurauia_tristyla"]
        self.assertEqual(taxon["include_gene_loss_denominator"], "false")
        self.assertEqual(taxon["include_species_tree"], "false")
        self.assertIn("Published Actinidiaceae", taxon["role"])
        self.assertIn("no S. tristyla assembled Angiosperms353 record", taxon["selection_note"])
        self.assertIn("Do not graft", taxon["selection_note"])
        self.assertEqual(
            self.brackets["Saurauia tristyla"]["current_status"],
            "prohibited_no_compatible_nuclear_asset_minimal_policy",
        )

    def test_one_species_one_leea_haplotype_contract_is_cross_registered(self) -> None:
        leea_rows = [row for row in self.taxa.values() if row["taxon_id"] == "leea_coccinea"]
        self.assertEqual(len(leea_rows), 1)
        leea = leea_rows[0]
        self.assertEqual(leea["biological_species"], "Leea coccinea")
        self.assertEqual(leea["include_gene_loss_denominator"], "false")
        self.assertIn("at most one QC-passing Lco1464 haplotype", leea["selection_note"])
        self.assertIn(
            "one declared Lco1464 representative haplotype",
            self.brackets["Leea coccinea"]["notes"],
        )
        self.assertIn("outside Vitoideae", self.brackets["Leea coccinea"]["notes"])

        catharanthus = self.taxa["catharanthus_roseus"]
        self.assertEqual(catharanthus["include_gene_loss_denominator"], "false")
        self.assertIn("GCA_024505715.1", catharanthus["selection_note"])
        self.assertEqual(
            self.brackets["Catharanthus roseus"]["current_status"],
            "exact_ncbi_matched_bundle_selected_not_a_rubiaceae_bracket",
        )

    def test_saurauia_rnaseq_registry_is_frozen_but_acquisition_is_prohibited(self) -> None:
        header, rows = read_tsv(CONFIG / "saurauia_nuclear_rescue_candidates.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["taxon_id"] for row in rows}, {"saurauia_tristyla"})
        self.assertEqual(
            {row["run_accession"] for row in rows},
            {"SRR11994221", "SRR28027655"},
        )
        self.assertEqual(
            {row["status"] for row in rows},
            {"prohibited_not_downloaded_minimal_outgroup_policy"},
        )
        self.assertTrue(all(row["layout"] == "paired" for row in rows))
        self.assertTrue(all(int(row["spots"]) > 20_000_000 for row in rows))
        self.assertTrue(all(int(row["bases"]) > 6_000_000_000 for row in rows))
        self.assertIn("planned_use", header)


if __name__ == "__main__":
    unittest.main()
