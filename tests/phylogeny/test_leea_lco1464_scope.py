"""Frozen same-individual haplotype scope for the Lco1464 outgroup."""

from __future__ import annotations

import csv
import re
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "phylogeny"
RECORD_COLUMNS = [
    "source_seqid",
    "gff_seqid",
    "canonical_seqid",
    "assembly_unit_id",
    "haplotype",
    "length_bp",
    "annotation_feature_rows",
    "gene_rows",
    "annotation_status",
    "sequence_class",
]
SUMMARY_COLUMNS = [
    "assembly_unit_id",
    "biological_species",
    "source_accession",
    "haplotype",
    "same_individual_group",
    "sequence_count",
    "assembly_bp",
    "n50_bp",
    "max_bp",
    "annotated_sequence_count",
    "annotated_sequence_bp",
    "annotation_feature_rows",
    "gene_loci",
    "primary_role",
    "chromosome_scope",
    "selection_status",
]
SEQID = re.compile(r"^Lco1464_v1\.0_(h[12])tg([0-9]{6})[lc]$")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


class LeeaLco1464ScopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record_header, cls.records = read_tsv(
            CONFIG / "leea_lco1464_records.tsv"
        )
        cls.summary_header, summaries = read_tsv(
            CONFIG / "leea_lco1464_haplotype_scope.tsv"
        )
        cls.summaries = {row["haplotype"]: row for row in summaries}

    def test_record_registry_is_exact_unique_and_does_not_fabricate_chromosomes(self) -> None:
        self.assertEqual(self.record_header, RECORD_COLUMNS)
        self.assertEqual(len(self.records), 140)
        source_ids = [row["source_seqid"] for row in self.records]
        self.assertEqual(len(source_ids), len(set(source_ids)))
        self.assertEqual(
            Counter(row["haplotype"] for row in self.records), {"h1": 63, "h2": 77}
        )
        for row in self.records:
            match = SEQID.fullmatch(row["source_seqid"])
            self.assertIsNotNone(match)
            assert match is not None
            self.assertEqual(match.group(1), row["haplotype"])
            self.assertEqual(row["canonical_seqid"], row["source_seqid"])
            self.assertEqual(
                row["sequence_class"],
                "publisher_hifiasm_contig_no_chromosome_assignment",
            )
            self.assertNotIn("Chr", row["canonical_seqid"])

    def test_gff_membership_and_integer_fields_are_fail_closed(self) -> None:
        for row in self.records:
            length = int(row["length_bp"])
            feature_rows = int(row["annotation_feature_rows"])
            genes = int(row["gene_rows"])
            self.assertGreater(length, 0)
            self.assertGreaterEqual(feature_rows, genes)
            self.assertGreaterEqual(genes, 0)
            if row["annotation_status"] == "annotated":
                self.assertEqual(row["gff_seqid"], row["source_seqid"])
                self.assertGreater(feature_rows, 0)
                self.assertGreater(genes, 0)
            elif row["annotation_status"] == "no_gff_features":
                self.assertEqual(row["gff_seqid"], "")
                self.assertEqual(feature_rows, 0)
                self.assertEqual(genes, 0)
            else:
                self.fail(f"Unexpected annotation status: {row['annotation_status']!r}")

    def test_summary_reconciles_exactly_to_record_registry(self) -> None:
        self.assertEqual(self.summary_header, SUMMARY_COLUMNS)
        self.assertEqual(set(self.summaries), {"h1", "h2"})
        expected = {
            "h1": {
                "count": 63,
                "bp": 554_723_631,
                "n50": 44_182_854,
                "max": 62_568_551,
                "annotated_count": 46,
                "annotated_bp": 553_353_975,
                "features": 311_977,
                "genes": 30_024,
            },
            "h2": {
                "count": 77,
                "bp": 550_183_858,
                "n50": 43_820_822,
                "max": 54_319_948,
                "annotated_count": 60,
                "annotated_bp": 548_536_384,
                "features": 277_521,
                "genes": 25_109,
            },
        }
        for haplotype, values in expected.items():
            records = [row for row in self.records if row["haplotype"] == haplotype]
            summary = self.summaries[haplotype]
            annotated = [
                row for row in records if row["annotation_status"] == "annotated"
            ]
            self.assertEqual(int(summary["sequence_count"]), len(records), values)
            self.assertEqual(
                int(summary["assembly_bp"]),
                sum(int(row["length_bp"]) for row in records),
            )
            self.assertEqual(int(summary["assembly_bp"]), values["bp"])
            self.assertEqual(int(summary["n50_bp"]), values["n50"])
            self.assertEqual(int(summary["max_bp"]), values["max"])
            self.assertEqual(int(summary["annotated_sequence_count"]), len(annotated))
            self.assertEqual(len(annotated), values["annotated_count"])
            self.assertEqual(
                int(summary["annotated_sequence_bp"]),
                sum(int(row["length_bp"]) for row in annotated),
            )
            self.assertEqual(int(summary["annotated_sequence_bp"]), values["annotated_bp"])
            self.assertEqual(
                int(summary["annotation_feature_rows"]),
                sum(int(row["annotation_feature_rows"]) for row in records),
            )
            self.assertEqual(int(summary["annotation_feature_rows"]), values["features"])
            self.assertEqual(
                int(summary["gene_loci"]),
                sum(int(row["gene_rows"]) for row in records),
            )
            self.assertEqual(int(summary["gene_loci"]), values["genes"])
            self.assertEqual(summary["same_individual_group"], "Lco1464_individual_1464")
            self.assertIn("no_publisher_chromosome_map", summary["chromosome_scope"])

        self.assertEqual(
            sum(int(row["length_bp"]) for row in self.records), 1_104_907_489
        )
        self.assertEqual(sum(int(row["gene_rows"]) for row in self.records), 55_133)

    def test_one_species_tip_and_full_haplotype_sensitivity_are_declared(self) -> None:
        h1 = self.summaries["h1"]
        h2 = self.summaries["h2"]
        self.assertEqual(h1["primary_role"], "provisional_primary_preference")
        self.assertEqual(h2["primary_role"], "required_diagnostic_and_representative_swap")
        self.assertIn("pending_per_haplotype", h1["selection_status"])
        self.assertIn("pending_per_haplotype", h2["selection_status"])

        _, taxa = read_tsv(CONFIG / "taxa.tsv")
        rows = [row for row in taxa if row["taxon_id"] == "leea_coccinea"]
        self.assertEqual(len(rows), 1)
        self.assertIn("at most one QC-passing Lco1464 haplotype", rows[0]["selection_note"])
        self.assertIn("cannot bracket crown Vitoideae", rows[0]["selection_note"])

        scope_doc = (ROOT / "docs" / "LEEA_LCO1464_SCOPE.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("The seventeen records without GFF features in each haplotype remain", scope_doc)
        self.assertIn("Do not select a haplotype because it yields a", scope_doc)
        self.assertIn("Retain both h1 and h2", scope_doc)
        self.assertIn("Never count h1 and h2 as two", scope_doc)

    def test_materializer_maps_reconcile_exactly_to_registry(self) -> None:
        expected_header = ["genome_seqid", "gff_seqid", "canonical_seqid"]
        for haplotype in ("h1", "h2"):
            header, rows = read_tsv(
                ROOT
                / "config"
                / "chromosome_maps"
                / f"leea_coccinea_lco1464_{haplotype}.full_haplotype.tsv"
            )
            self.assertEqual(header, expected_header)
            registry_rows = [
                row for row in self.records if row["haplotype"] == haplotype
            ]
            self.assertEqual(len(rows), len(registry_rows))
            self.assertEqual(
                rows,
                [
                    {
                        "genome_seqid": row["source_seqid"],
                        "gff_seqid": row["gff_seqid"],
                        "canonical_seqid": row["canonical_seqid"],
                    }
                    for row in registry_rows
                ],
            )
            self.assertEqual(sum(row["gff_seqid"] == "" for row in rows), 17)


if __name__ == "__main__":
    unittest.main()
