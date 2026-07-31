"""Static contract tests for the frozen *Actinidia rufa* identity audit."""

from __future__ import annotations

import csv
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "config" / "rufa_bundle_scope_identity.tsv"
PAIRWISE = ROOT / "config" / "rufa_bundle_pairwise_identity.tsv"
MAP = (
    ROOT
    / "config"
    / "chromosome_maps"
    / "act_rufa_actinidiabase_v1.publisher_scope.tsv"
)
ARU_MAP = (
    ROOT / "config" / "chromosome_maps" / "act_rufa_aru_r1.publisher_scope.tsv"
)
FUCHU_MAP = (
    ROOT / "config" / "chromosome_maps" / "act_rufa_fuchu.publisher_scope.tsv"
)
EXCLUDED = (
    ROOT
    / "config"
    / "chromosome_maps"
    / "act_rufa_actinidiabase_v1.excluded_records.tsv"
)
ASSET_CHECKSUMS = ROOT / "config" / "rufa_bundle_asset_checksums.tsv"
ASSEMBLIES = ROOT / "config" / "assemblies.tsv"
DOWNLOADS = ROOT / "config" / "downloads.tsv"
CATALOG = ROOT / "config" / "candidate_assembly_catalog.tsv"

EXPECTED_UNITS = {
    "act_rufa_aru_r1",
    "act_rufa_fuchu",
    "act_rufa_actinidiabase_v1",
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class RufaScopeIdentityTest(unittest.TestCase):
    def test_summary_scope_closes_exactly(self) -> None:
        observed = {row["assembly_unit_id"]: row for row in rows(SUMMARY)}
        self.assertEqual(set(observed), EXPECTED_UNITS)

        expected = {
            "act_rufa_aru_r1": (29, 620324227, 29, 620324227, 0, 0, 29, 658027),
            "act_rufa_fuchu": (
                501,
                647239177,
                29,
                509119555,
                472,
                138119622,
                501,
                299762,
            ),
            "act_rufa_actinidiabase_v1": (
                38,
                615891845,
                29,
                613679866,
                9,
                2211979,
                37,
                533866,
            ),
        }
        numeric_fields = (
            "genome_records",
            "genome_bp",
            "publisher_chromosome_records",
            "publisher_chromosome_bp",
            "extra_genome_records",
            "extra_genome_bp",
            "gff_seqids",
            "gff_feature_rows",
        )
        for unit, row in observed.items():
            self.assertEqual(
                tuple(int(row[field]) for field in numeric_fields), expected[unit]
            )
            self.assertEqual(
                int(row["genome_records"]),
                int(row["publisher_chromosome_records"])
                + int(row["extra_genome_records"]),
            )
            self.assertEqual(
                int(row["genome_bp"]),
                int(row["publisher_chromosome_bp"])
                + int(row["extra_genome_bp"]),
            )
            self.assertEqual(
                int(row["gff_seqids"]),
                int(row["chromosome_gff_seqids"])
                + int(row["extra_gff_seqids"]),
            )
            self.assertEqual(
                int(row["gff_feature_rows"]),
                int(row["chromosome_gff_feature_rows"])
                + int(row["extra_gff_feature_rows"]),
            )
            for field in (
                "genome_gzip_sha256",
                "gff_gzip_sha256",
                "whole_sequence_multiset_sha256",
                "chromosome_sequence_multiset_sha256",
            ):
                self.assertRegex(row[field], HEX64)

    def test_all_pairs_are_sequence_distinct_without_claiming_biological_identity(self) -> None:
        observed = rows(PAIRWISE)
        self.assertEqual(len(observed), 3)
        pairs = {
            frozenset((row["assembly_unit_id_a"], row["assembly_unit_id_b"]))
            for row in observed
        }
        self.assertEqual(
            pairs,
            {
                frozenset((a, b))
                for index, a in enumerate(sorted(EXPECTED_UNITS))
                for b in sorted(EXPECTED_UNITS)[index + 1 :]
            },
        )
        for row in observed:
            self.assertEqual(row["exact_sequence_record_matches"], "0")
            self.assertEqual(row["exact_matching_bp"], "0")
            self.assertEqual(
                row["sequence_identity_conclusion"],
                "sequence_distinct_not_exact_mirrors",
            )
        actinidiabase_pairs = [
            row
            for row in observed
            if "act_rufa_actinidiabase_v1"
            in {row["assembly_unit_id_a"], row["assembly_unit_id_b"]}
        ]
        self.assertTrue(
            all(
                row["biological_identity_conclusion"]
                == "unresolved_ActinidiaBase_biological_accession"
                for row in actinidiabase_pairs
            )
        )

    def test_actinidiabase_paper_contig_metric_is_not_release_record_count(self) -> None:
        observed = {row["assembly_unit_id"]: row for row in rows(SUMMARY)}
        actinidiabase = observed["act_rufa_actinidiabase_v1"]
        self.assertEqual(actinidiabase["publication_reported_contig_count"], "100")
        self.assertEqual(actinidiabase["publication_reported_assembly_bp"], "615891845")
        self.assertEqual(actinidiabase["genome_records"], "38")
        self.assertEqual(actinidiabase["genome_bp"], "615891845")
        self.assertIn("supplementary_table_S1", actinidiabase["publication_metric_source"])
        for unit in ("act_rufa_aru_r1", "act_rufa_fuchu"):
            self.assertEqual(
                observed[unit]["publication_reported_contig_count"], "not_audited"
            )

    def test_actinidiabase_scope_map_is_explicit_and_complete(self) -> None:
        observed = rows(MAP)
        self.assertEqual(len(observed), 29)
        self.assertEqual(
            list(observed[0]), ["genome_seqid", "gff_seqid", "canonical_seqid"]
        )
        self.assertEqual(
            [row["genome_seqid"] for row in observed],
            [f"Chr{index}" for index in range(1, 30)],
        )
        self.assertEqual(
            [row["gff_seqid"] for row in observed],
            [f"Chr{index}" for index in range(1, 30)],
        )
        self.assertEqual(
            [row["canonical_seqid"] for row in observed],
            [f"PubChr{index:02d}" for index in range(1, 30)],
        )

        excluded = rows(EXCLUDED)
        self.assertEqual(len(excluded), 9)
        self.assertEqual(sum(int(row["length_bp"]) for row in excluded), 2211979)
        self.assertEqual(sum(int(row["gff_feature_rows"]) for row in excluded), 1642)
        no_features = [row for row in excluded if row["gff_feature_rows"] == "0"]
        self.assertEqual([row["genome_seqid"] for row in no_features], ["Contig01298"])
        self.assertTrue(all(HEX64.fullmatch(row["sequence_sha256"]) for row in excluded))

    def test_named_aru_and_fuchu_scope_maps_match_exact_namespaces(self) -> None:
        aru = rows(ARU_MAP)
        fuchu = rows(FUCHU_MAP)
        self.assertEqual(len(aru), 29)
        self.assertEqual(len(fuchu), 29)
        for index, row in enumerate(aru, start=1):
            self.assertEqual(row["genome_seqid"], f"BRYG010000{index:02d}.1")
            self.assertEqual(row["gff_seqid"], f"ARU1.0ch{index:02d}")
            self.assertEqual(row["canonical_seqid"], f"PubChr{index:02d}")
        for index, row in enumerate(fuchu, start=1):
            self.assertEqual(row["genome_seqid"], f"BJWL010000{index:02d}.1")
            self.assertEqual(row["gff_seqid"], f"BJWL010000{index:02d}.1")
            self.assertEqual(row["canonical_seqid"], f"PubChr{index:02d}")

    def test_executable_config_and_download_hashes_match_audit(self) -> None:
        summary = {row["assembly_unit_id"]: row for row in rows(SUMMARY)}
        assemblies = {
            row["assembly_unit_id"]: row
            for row in rows(ASSEMBLIES)
            if row["assembly_unit_id"] in EXPECTED_UNITS
        }
        self.assertEqual(set(assemblies), EXPECTED_UNITS)
        for unit, row in assemblies.items():
            self.assertEqual(
                row["expected_genome_sha256"], summary[unit]["genome_gzip_sha256"]
            )
            self.assertEqual(
                row["expected_annotation_sha256"], summary[unit]["gff_gzip_sha256"]
            )
            self.assertEqual(row["partition_rule"], "explicit_publisher_scope_map")

        rufa_downloads = [
            row for row in rows(DOWNLOADS) if row["assembly_unit_id"] in EXPECTED_UNITS
        ]
        self.assertEqual(len(rufa_downloads), 12)
        audit_assets = {row["asset_id"]: row for row in rows(ASSET_CHECKSUMS)}
        self.assertEqual(set(audit_assets), {row["asset_id"] for row in rufa_downloads})
        for row in rufa_downloads:
            audit = audit_assets[row["asset_id"]]
            self.assertEqual(audit["assembly_unit_id"], row["assembly_unit_id"])
            self.assertEqual(audit["asset_type"], row["asset_type"])
            self.assertEqual(audit["bytes"], row["expected_bytes"])
            self.assertRegex(audit["sha256"], HEX64)
            self.assertEqual(
                audit["publisher_md5"], row["md5"] if row["md5"] else "none"
            )
            if row["asset_type"] == "genome":
                self.assertEqual(
                    audit["sha256"],
                    summary[row["assembly_unit_id"]]["genome_gzip_sha256"],
                )
            if row["asset_type"] == "gff":
                self.assertEqual(
                    audit["sha256"],
                    summary[row["assembly_unit_id"]]["gff_gzip_sha256"],
                )

    def test_catalog_does_not_call_actinidiabase_an_aru_or_fuchu_mirror(self) -> None:
        catalog = {
            row["catalog_id"]: row for row in rows(CATALOG)
        }["act_rufa_actinidiabase_2024"]
        self.assertEqual(
            catalog["assembly_scope"], "38_records_29_pseudochromosomes_plus_9_contigs"
        )
        self.assertEqual(
            catalog["biological_independence"],
            "sequence_distinct_candidate_biological_accession_unresolved",
        )
        self.assertIn("zero exact normalized sequence records", catalog["decision_note"])
        self.assertIn("do not count", catalog["decision_note"])


if __name__ == "__main__":
    unittest.main()
