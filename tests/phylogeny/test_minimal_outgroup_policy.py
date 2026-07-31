"""Fail-closed contracts for the approved minimal-outgroup and dating design."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "phylogeny"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


class MinimalOutgroupPolicyTest(unittest.TestCase):
    def test_primary_membership_is_bounded_and_uses_no_new_taxa(self) -> None:
        expected_header = [
            "outgroup_set_id",
            "entry_id",
            "eligible_taxon_ids",
            "min_tips",
            "max_tips",
            "relation_to_primary",
            "selection_rule",
            "phylogenetic_role",
            "calibration_effect",
            "download_policy",
            "status",
            "evidence_source",
        ]
        header, rows = read_tsv(CONFIG / "minimal_outgroup_design.tsv")
        self.assertEqual(header, expected_header)
        self.assertEqual(
            len({(row["outgroup_set_id"], row["entry_id"]) for row in rows}),
            len(rows),
        )

        primary = {
            row["entry_id"]: row
            for row in rows
            if row["outgroup_set_id"] == "primary_minimal"
        }
        self.assertEqual(
            set(primary),
            {
                "one_rhododendron",
                "coffea_arabica",
                "vitis_vinifera",
                "leea_coccinea",
                "catharanthus_roseus",
            },
        )
        self.assertEqual(
            (primary["one_rhododendron"]["min_tips"], primary["one_rhododendron"]["max_tips"]),
            ("1", "1"),
        )
        self.assertEqual(
            (primary["coffea_arabica"]["min_tips"], primary["coffea_arabica"]["max_tips"]),
            ("1", "1"),
        )
        self.assertEqual(
            (primary["vitis_vinifera"]["min_tips"], primary["vitis_vinifera"]["max_tips"]),
            ("1", "1"),
        )
        self.assertEqual(
            (primary["leea_coccinea"]["min_tips"], primary["leea_coccinea"]["max_tips"]),
            ("0", "1"),
        )
        self.assertEqual(
            (primary["catharanthus_roseus"]["min_tips"], primary["catharanthus_roseus"]["max_tips"]),
            ("0", "0"),
        )
        self.assertTrue(
            all(
                row["download_policy"]
                in {
                    "existing_assets_only_no_new_taxa",
                    "already_acquired_no_new_taxa",
                    "no_further_acquisition",
                }
                for row in rows
            )
        )

    def test_leea_catharanthus_and_saurauia_roles_are_not_interchangeable(self) -> None:
        _, rows = read_tsv(CONFIG / "minimal_outgroup_design.tsv")
        primary = {
            row["entry_id"]: row
            for row in rows
            if row["outgroup_set_id"] == "primary_minimal"
        }
        self.assertIn("outside Vitoideae", primary["leea_coccinea"]["calibration_effect"])
        self.assertIn("cannot bracket", primary["catharanthus_roseus"]["calibration_effect"])

        cath_swap = [
            row
            for row in rows
            if row["outgroup_set_id"] == "catharanthus_for_coffea_root_swap"
        ]
        self.assertEqual(len(cath_swap), 1)
        self.assertEqual(cath_swap[0]["relation_to_primary"], "replaces_coffea_arabica")

        saurauia = [row for row in rows if row["entry_id"] == "saurauia_tristyla"]
        self.assertEqual(len(saurauia), 1)
        self.assertEqual((saurauia[0]["min_tips"], saurauia[0]["max_tips"]), ("0", "0"))
        self.assertEqual(saurauia[0]["status"], "prohibited_no_compatible_nuclear_genome_or_proteome")

    def test_corrected_fossil_nodes_remain_disabled(self) -> None:
        _, rows = read_tsv(CONFIG / "calibrations.tsv")
        by_id = {row["calibration_id"]: row for row in rows}
        self.assertNotIn("vitaceae_crown_indovitis", by_id)

        indovitis = by_id["vitoideae_crown_indovitis"]
        self.assertEqual(indovitis["calibrated_clade"], "Vitoideae")
        self.assertEqual(indovitis["primary_source"], "https://doi.org/10.3732/ajb.1300008")
        self.assertIn("Leea is the Leeoideae outgroup", indovitis["required_taxon_sampling"])
        self.assertIn("two sampled Vitoideae", indovitis["required_taxon_sampling"])
        self.assertEqual(indovitis["status"], "disabled_missing_second_vitoideae_lineage")

        rubiaceae = by_id["rubiaceae_internal_candidate"]
        self.assertEqual(rubiaceae["minimum_ma"], "")
        self.assertEqual(rubiaceae["maximum_ma"], "")
        self.assertIn("Catharanthus is Apocynaceae", rubiaceae["required_taxon_sampling"])
        self.assertTrue(rubiaceae["status"].startswith("disabled_"))
        self.assertFalse(any(row["status"] == "active" for row in rows))

    def test_unacquired_candidates_are_disabled_in_taxon_registry(self) -> None:
        _, rows = read_tsv(CONFIG / "taxa.tsv")
        by_id = {row["taxon_id"]: row for row in rows}
        for taxon_id in {
            "coffea_canephora",
            "rubia_or_cinchona",
            "non_rhododendron_ericaceae",
            "actinidiaceae_sister_outgroup",
            "ampelocissus_or_nothocissus",
        }:
            self.assertEqual(by_id[taxon_id]["include_species_tree"], "false")
            self.assertEqual(
                by_id[taxon_id]["current_asset_status"],
                "not_acquired_minimal_outgroup_policy",
            )

        self.assertEqual(by_id["saurauia_tristyla"]["include_species_tree"], "false")
        _, rescue_rows = read_tsv(CONFIG / "saurauia_nuclear_rescue_candidates.tsv")
        self.assertEqual(
            {row["status"] for row in rescue_rows},
            {"prohibited_not_downloaded_minimal_outgroup_policy"},
        )

    def test_dating_designs_separate_primary_secondary_and_legacy(self) -> None:
        expected_ids = {
            "revised_primary_fossil",
            "secondary_timetree_sensitivity",
            "legacy_timetree_chronos_reproduction",
        }
        _, rows = read_tsv(CONFIG / "dating_designs.tsv")
        by_id = {row["dating_design_id"]: row for row in rows}
        self.assertEqual(set(by_id), expected_ids)
        self.assertIn("blocked_no_active", by_id["revised_primary_fossil"]["status"])
        self.assertIn("undocumented point ages are prohibited", by_id["secondary_timetree_sensitivity"]["age_encoding"])
        self.assertTrue(by_id["secondary_timetree_sensitivity"]["status"].startswith("active_"))
        self.assertIn("TimeTree secondary-calibrated", by_id["secondary_timetree_sensitivity"]["publication_claim"])
        self.assertIn("time_min equals time_max", by_id["legacy_timetree_chronos_reproduction"]["age_encoding"])
        self.assertIn("9_tip", by_id["legacy_timetree_chronos_reproduction"]["status"])

        secondary_header, secondary_rows = read_tsv(CONFIG / "secondary_timetree_constraints.tsv")
        self.assertEqual(len(secondary_rows), 4)
        self.assertTrue(all(row["status"].startswith("active_") for row in secondary_rows))
        self.assertEqual(
            {row["constraint_id"] for row in secondary_rows},
            {
                "timetree_vitis_coffea_root",
                "timetree_coffea_rhododendron",
                "timetree_rhododendron_actinidia",
                "timetree_actinidia_crown",
            },
        )
        for required in {
            "timetree_version",
            "retrieved_utc",
            "query_url",
            "contributing_studies",
            "raw_artifact_sha256",
            "transformation_note",
            "status",
        }:
            self.assertIn(required, secondary_header)

    def test_exact_legacy_points_and_broken_provenance_are_frozen(self) -> None:
        _, rows = read_tsv(CONFIG / "legacy_chronos_calibrations.tsv")
        self.assertEqual(len(rows), 5)
        self.assertEqual(
            {row["script_literal_ma"] for row in rows},
            {"113.97", "41.86", "5.91", "2.29", "25.33"},
        )
        self.assertEqual(
            {row["timetree_artifact_ma"] for row in rows},
            {"113.97482", "41.86381", "5.91334", "2.29026", "25.33251"},
        )
        self.assertEqual(
            {row["script_sha256"] for row in rows},
            {"290f926fb3f7e6bc3822278df3d0809799de6c4da4c93b11e90437ab71b08929"},
        )
        self.assertEqual(
            {row["timetree_sha256"] for row in rows},
            {"ac63bef64c87a9b18f42e01662060dc58352009564fedb06424f8c31154dc9d3"},
        )
        self.assertEqual({row["production_use"] for row in rows}, {"prohibited_in_revised_production"})
        self.assertTrue(all("time_min equal to time_max" in row["source_state"] for row in rows))

        _, artifacts = read_tsv(CONFIG / "legacy_chronos_artifacts.tsv")
        by_id = {row["artifact_id"]: row for row in artifacts}
        self.assertEqual(set(by_id), {"legacy_calibrate_script", "legacy_timetree_points", "current_script_input", "archived_calibrated_output"})
        self.assertEqual(by_id["current_script_input"]["tip_count"], "9")
        self.assertEqual(by_id["archived_calibrated_output"]["tip_count"], "17")
        self.assertIn("fail_tip_count_mismatch", by_id["current_script_input"]["audit_status"])
        self.assertEqual({row["production_use"] for row in artifacts}, {"prohibited"})

    def test_policy_document_explains_why_and_how(self) -> None:
        text = (ROOT / "docs" / "OUTGROUPS_AND_CALIBRATIONS.md").read_text(encoding="utf-8")
        for phrase in (
            "No additional outgroup taxon is downloaded",
            "one-for-one replacement for *Coffea*",
            "cannot define crown Vitoideae",
            "cannot define a node internal to Rubiaceae",
            "Three non-interchangeable dating designs",
            "9-tip tree",
            "17 tips",
            "no undocumented fixed point",
        ):
            self.assertIn(phrase, text)


if __name__ == "__main__":
    unittest.main()
