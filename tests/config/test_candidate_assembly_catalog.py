"""Fail-closed validation for the public candidate assembly inventory."""

from __future__ import annotations

import csv
import re
import unittest
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "config" / "candidate_assembly_catalog.tsv"

EXPECTED_COLUMNS = [
    "catalog_id",
    "biological_species",
    "biological_accession",
    "assembly_unit_label",
    "year",
    "assembly_accession",
    "alternate_accessions",
    "accession_version_status",
    "bioproject",
    "ploidy_or_representation",
    "assembly_scope",
    "public_genome",
    "public_gff",
    "public_cds",
    "public_protein",
    "classification",
    "deduplication_group",
    "biological_independence",
    "repository_url",
    "publication_doi",
    "decision_note",
]


def read_catalog() -> list[dict[str, str]]:
    with CATALOG.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise AssertionError(
                f"Unexpected candidate catalog columns: {reader.fieldnames!r}"
            )
        rows = list(reader)
    if not rows:
        raise AssertionError("Candidate assembly catalog is empty")
    return rows


class CandidateAssemblyCatalogTest(unittest.TestCase):
    def setUp(self) -> None:
        self.rows = read_catalog()
        self.by_id = {row["catalog_id"]: row for row in self.rows}

    def test_schema_values_and_public_safety(self) -> None:
        self.assertEqual(len(self.by_id), len(self.rows), "catalog_id is not unique")
        self.assertEqual(
            {row["classification"] for row in self.rows},
            {
                "primary_replacement_candidate",
                "sensitivity",
                "unavailable",
                "excluded",
            },
        )
        for row in self.rows:
            self.assertTrue(all(row[column].strip() for column in EXPECTED_COLUMNS))
            self.assertRegex(row["year"], r"^\d{4}$")
            for column in ("public_genome", "public_gff", "public_cds", "public_protein"):
                self.assertIn(row[column], {"yes", "no"})
            repository_url = row["repository_url"]
            self.assertTrue(
                repository_url == "not_available" or repository_url.startswith("https://")
            )
            joined = "\t".join(row.values())
            self.assertTrue(all(not value.startswith("/") for value in row.values()))
            credential_terms = ("to" + "ken", "sec" + "ret", "pass" + "word")
            self.assertFalse(
                any(re.search(rf"{term}=", joined, re.I) for term in credential_terms)
            )

    def test_classification_is_fail_closed(self) -> None:
        for row in self.rows:
            runnable = all(
                row[column] == "yes"
                for column in ("public_genome", "public_gff", "public_protein")
            )
            if row["classification"] == "primary_replacement_candidate":
                self.assertTrue(runnable, row["catalog_id"])
            if row["classification"] == "unavailable":
                self.assertFalse(runnable, row["catalog_id"])

        dna_only = [
            row
            for row in self.rows
            if row["public_genome"] == "yes"
            and row["public_gff"] == "no"
            and row["public_protein"] == "no"
        ]
        self.assertGreaterEqual(len(dna_only), 20)
        self.assertTrue(all(row["classification"] == "unavailable" for row in dna_only))

    def test_paired_haplotypes_are_one_biological_accession(self) -> None:
        paired = defaultdict(list)
        for row in self.rows:
            if row["biological_independence"] in {
                "paired_haplotype_same_accession",
                "hybrid_haplotype_same_accession",
            }:
                paired[row["deduplication_group"]].append(row)

        self.assertGreaterEqual(len(paired), 20)
        for group, rows in paired.items():
            self.assertEqual(len(rows), 2, group)
            self.assertEqual(len({row["biological_accession"] for row in rows}), 1, group)
            self.assertTrue(
                all("independent species" in row["decision_note"] for row in rows),
                group,
            )

    def test_required_biological_and_availability_warnings(self) -> None:
        hap1 = self.by_id["act_eriantha_hap1_2026"]
        hap2 = self.by_id["act_eriantha_hap2_2026"]
        self.assertEqual(hap1["biological_accession"], hap2["biological_accession"])
        self.assertEqual(hap1["deduplication_group"], hap2["deduplication_group"])
        self.assertTrue(
            all(
                row["assembly_scope"] == "29_anchored_pseudochromosomes_only"
                and "chromosome-anchored subset" in row["decision_note"]
                for row in (hap1, hap2)
            )
        )

        hh01 = [row for row in self.rows if row["biological_accession"] == "HH01"]
        self.assertEqual(len(hh01), 2)
        self.assertTrue(all(row["classification"] == "excluded" for row in hh01))
        self.assertTrue(all("hybrid" in row["decision_note"].lower() for row in hh01))

        zhejiangensis = [
            row for row in self.rows if row["biological_accession"] == "legacy_A_zhejiangensis"
        ]
        self.assertEqual(len(zhejiangensis), 2)
        self.assertEqual(
            {row["assembly_unit_label"] for row in zhejiangensis},
            {"parental_haplome_A", "parental_haplome_B"},
        )
        self.assertTrue(all(row["classification"] == "excluded" for row in zhejiangensis))
        self.assertTrue(
            all("not independent species" in row["decision_note"] for row in zhejiangensis)
        )

        mt = self.by_id["act_rufa_mt570001_2026"]
        self.assertEqual(mt["classification"], "unavailable")
        self.assertEqual(mt["assembly_accession"], "not_available")
        self.assertEqual(
            Counter(mt[column] for column in ("public_genome", "public_gff", "public_cds", "public_protein")),
            Counter({"no": 4}),
        )

        guimi = self.by_id["act_chinensis_guimi2_2025"]
        self.assertEqual(
            set(guimi["alternate_accessions"].split(";")),
            {"GCA_051903715.1", "GCA_051167385.1"},
        )
        self.assertIn("must be counted once", guimi["decision_note"])
        self.assertEqual(guimi["classification"], "unavailable")
        self.assertEqual(
            tuple(
                guimi[column]
                for column in ("public_genome", "public_gff", "public_cds", "public_protein")
            ),
            ("yes", "no", "no", "no"),
        )
        self.assertIn("no matched structural annotation", guimi["decision_note"])

        for catalog_id in ("act_chinensis_kuimi_2026", "act_chinensis_acm4_2026"):
            row = self.by_id[catalog_id]
            self.assertEqual(row["ploidy_or_representation"], "tetraploid_4x")
            self.assertEqual(row["assembly_scope"], "chromosome_level_116_chromosomes")
            self.assertEqual(row["classification"], "sensitivity")
            self.assertIn("ploidy-mismatched sensitivity", row["decision_note"])


if __name__ == "__main__":
    unittest.main()
