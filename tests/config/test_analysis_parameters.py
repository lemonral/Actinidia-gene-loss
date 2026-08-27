"""Regression checks for frozen project-wide analysis parameters."""

from __future__ import annotations

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class AnalysisParameterTest(unittest.TestCase):
    def test_worker_and_jcvi_policy_is_explicit(self) -> None:
        data = tomllib.loads(
            (ROOT / "config" / "analysis_parameters.toml").read_text(encoding="utf-8")
        )
        self.assertEqual(data["maximum_project_workers"], 15)
        jcvi = data["assembly_qc"]["jcvi"]
        self.assertEqual(jcvi["primary_gate_minimum_fraction"], 0.50)
        self.assertTrue(jcvi["require_reference_and_query_coverage_report"])
        self.assertEqual(jcvi["aligner"], "LAST")
        self.assertEqual(jcvi["database_type"], "protein")
        self.assertEqual(jcvi["cscore"], 0.7)
        self.assertEqual(jcvi["tandem_nmax"], 10)
        self.assertEqual(jcvi["maximum_gene_distance"], 20)
        self.assertEqual(jcvi["minimum_anchor_block_size"], 4)
        self.assertEqual(jcvi["coverage_anchor_source"], "raw JCVI anchors")
        self.assertEqual(data["gene_loss"]["shared_definition"],
                         "decayed or deleted in all 23 included genomes")
        self.assertFalse(data["gene_loss"]["uncertain_is_positive"])
        self.assertEqual(
            data["gene_loss"]["positive_classes"],
            ["decayed", "deleted"],
        )
        self.assertEqual(
            data["gene_loss"]["strict_sensitivity_positive_classes"],
            ["pseudogenized", "deleted"],
        )
        self.assertEqual(
            data["gene_loss"]["comparative_analysis_set"],
            "shared_and_non_shared",
        )
        self.assertFalse(data["gene_loss"]["uncertain_enters_resolved_denominator"])
        self.assertEqual(
            data["gene_loss"]["historical_reproduction_positive_classes"],
            ["decayed", "deleted"],
        )
        self.assertEqual(data["gene_loss"]["tblastx"]["maximum_workers"], 14)
        miniprot = data["gene_loss"]["miniprot"]
        self.assertEqual(miniprot["maximum_workers"], 14)
        self.assertEqual(miniprot["version"], "0.18-r281")
        self.assertEqual(miniprot["strict_disruption_tags"], ["fs", "st"])
        self.assertEqual(miniprot["strict_disruption_query_coverage_minimum"], 0.80)
        self.assertEqual(miniprot["strict_disruption_exact_identity_minimum"], 0.70)
        self.assertEqual(miniprot["strict_disruption_alignment_score_minimum"], 100)
        self.assertEqual(
            data["spatial"]["primary_positive_classes"],
            ["decayed"],
        )
        self.assertEqual(
            data["spatial"]["resolved_denominator_classes"],
            [],
        )
        self.assertIn(
            "annotated target genes",
            data["spatial"]["primary_denominator"],
        )
        self.assertEqual(
            data["spatial"]["deleted_position_role"],
            "excluded because no observed residual-sequence locus is available",
        )
        self.assertEqual(data["nlr"]["maximum_workers"], 8)
        self.assertTrue(data["nlr"]["shared_positive_complete_genes_excluded"])
        self.assertFalse(data["nlr"]["uncertain_callable_comparisons_remain_in_denominator"])
        homology = data["chromosome_homology"]
        self.assertIn("not the chromosome-naming authority", homology["production_role"])
        self.assertEqual(homology["matrix_schema_version"], "1.0.0")
        self.assertEqual(homology["matrix_provenance_schema_version"], "1.0.0")
        self.assertRegex(homology["reference_asset_registry_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(
            homology["reference_chromosome_map_registry_sha256"],
            r"^[0-9a-f]{64}$",
        )
        self.assertEqual(homology["expected_query_chromosomes"], 29)
        self.assertEqual(homology["expected_reference_chromosomes"], 29)
        self.assertEqual(homology["minimap2_version"], "2.28-r1209")
        self.assertEqual(
            homology["minimap2_command_template"],
            [
                "minimap2",
                "-x",
                "asm5",
                "--secondary=no",
                "-c",
                "--cs=long",
                "{reference_fasta}",
                "{query_fasta}",
            ],
        )
        self.assertEqual(
            homology["nucleotide_score_formula"],
            "harmonic_mean(query_coverage,reference_coverage)*(1-weighted_divergence)",
        )
        self.assertEqual(
            homology["reciprocal_nucleotide_coverage_formula"],
            "min(query_coverage,reference_coverage)",
        )
        self.assertEqual(
            homology["jcvi_score_formula"],
            "harmonic_mean(query_gene_coverage,reference_gene_coverage)",
        )
        self.assertTrue(homology["minimap2_primary_only"])
        self.assertEqual(homology["minimum_mapq"], 20)
        self.assertEqual(homology["minimum_alignment_block_bp"], 10000)
        self.assertEqual(homology["maximum_de"], 0.15)
        self.assertEqual(homology["minimum_top_second_ratio"], 1.5)
        self.assertEqual(homology["minimum_normalized_score_margin"], 0.10)
        self.assertEqual(
            homology["minimum_assigned_reciprocal_nucleotide_coverage"], 0.05
        )
        self.assertEqual(homology["minimum_unique_anchor_pairs"], 30)
        self.assertEqual(homology["minimum_reciprocal_gene_coverage"], 0.05)
        self.assertEqual(homology["minimum_assigned_jcvi_score"], 0.05)
        self.assertEqual(homology["arithmetic_tolerance"], 1e-9)
        self.assertTrue(homology["require_row_and_column_reciprocal_best"])
        self.assertTrue(homology["require_nucleotide_and_jcvi_assignment_agreement"])
        self.assertTrue(homology["require_hy4a_hy4p_label_agreement"])
        self.assertFalse(homology["reverse_complement_allowed"])
        orientation = data["chromosome_orientation_harmonization"]
        self.assertIn("superseded", orientation["production_status"])
        self.assertEqual(
            orientation["production_direction_policy"],
            "preserve publisher direction; relabel only",
        )
        self.assertEqual(
            orientation["primary_orientation_reference"],
            "act_chinensis_hongyang_v4_hy4a",
        )
        self.assertEqual(orientation["minimum_dominant_orientation_fraction"], 0.80)
        self.assertEqual(orientation["minimum_oriented_matching_bases"], 1_000_000)
        self.assertEqual(orientation["negative_action"], "reverse_complement")
        self.assertTrue(orientation["transform_gff_coordinates_and_strands"])
        self.assertTrue(
            orientation["require_cds_and_protein_sequence_identity_after_transform"]
        )
        self.assertFalse(data["release"]["push_automatically"])

    def test_translated_search_thresholds_are_frozen(self) -> None:
        data = tomllib.loads(
            (ROOT / "config" / "analysis_parameters.toml").read_text(encoding="utf-8")
        )
        policy = data["gene_loss"]["tblastx"]
        self.assertEqual(policy["postfilter_percent_identity_minimum"], 50.0)
        self.assertEqual(policy["postfilter_bitscore_minimum"], 50.0)
        self.assertEqual(policy["postfilter_evalue"], 1e-5)
        self.assertEqual(policy["postfilter_alignment_length_minimum"], 0)


if __name__ == "__main__":
    unittest.main()
