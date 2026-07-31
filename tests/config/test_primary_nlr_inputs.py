"""Contract tests for the primary NLR input cohort."""

from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PrimaryNlrInputConfigTest(unittest.TestCase):
    def test_target_set_exactly_matches_primary_loss_units(self) -> None:
        nlr = read_tsv(ROOT / "config" / "primary_nlr_input_sources.tsv")
        primary = read_tsv(ROOT / "config" / "primary_species_loss_aggregation.tsv")
        references = [row for row in nlr if row["analysis_role"] == "reference_callable"]
        targets = [row for row in nlr if row["analysis_role"] == "target_repertoire"]
        self.assertEqual(len(references), 1)
        self.assertEqual(len(targets), 23)
        self.assertEqual(
            {row["sample_id"] for row in targets},
            {row["assembly_unit_id"] for row in primary if row["include"] == "true"},
        )
        self.assertEqual({row["input_scope"] for row in targets}, {"whole_genome"})
        self.assertEqual({row["expected_fasta_records"] for row in targets}, {"29"})

    def test_output_names_and_sources_are_unique_relative_paths(self) -> None:
        rows = read_tsv(ROOT / "config" / "primary_nlr_input_sources.tsv")
        outputs = [row["output_basename"] for row in rows]
        sources = [row["source_fasta"] for row in rows]
        self.assertEqual(len(outputs), len(set(outputs)))
        self.assertEqual(len(sources), len(set(sources)))
        for value in sources:
            path = Path(value)
            self.assertFalse(path.is_absolute())
            self.assertNotIn("..", path.parts)

    def test_summary_metadata_matches_target_identity_and_preserves_unit_suffixes(self) -> None:
        nlr = read_tsv(ROOT / "config" / "primary_nlr_input_sources.tsv")
        metadata = read_tsv(ROOT / "config" / "primary_nlr_unit_metadata.tsv")
        targets = {
            row["sample_id"]: row for row in nlr if row["analysis_role"] == "target_repertoire"
        }
        self.assertEqual(len(metadata), 23)
        self.assertEqual(set(targets), {row["assembly_unit_id"] for row in metadata})
        self.assertEqual({row["include"] for row in metadata}, {"true"})
        self.assertEqual(
            {row["analysis_cohort"] for row in metadata}, {"primary_23_units_nonshared_v1"}
        )
        self.assertEqual(
            {row["assembly_scope"] for row in metadata}, {"primary_29_chromosome_analysis_unit"}
        )
        for row in metadata:
            self.assertEqual(row["biological_species"], targets[row["assembly_unit_id"]]["species"])
        suffixes = {
            row["assembly_unit_id"]: row["haplotype_or_subgenome"] for row in metadata
        }
        self.assertEqual([suffixes[f"act_arguta_{letter.lower()}_legacy"] for letter in "ABCD"], list("ABCD"))
        self.assertEqual(
            [suffixes[f"act_deliciosa_adm_2026_{letter}"] for letter in "ABCDEF"], list("ABCDEF")
        )
        self.assertEqual(suffixes["act_eriantha_hap1_2026"], "HAP1")
        self.assertEqual(suffixes["act_eriantha_hap2_2026"], "HAP2")


if __name__ == "__main__":
    unittest.main()
