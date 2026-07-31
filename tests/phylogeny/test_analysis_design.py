"""Regression checks for the phylogeny and gene-family design contracts."""

from __future__ import annotations

import csv
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 compatibility on the analysis server.
    import tomli as tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "phylogeny"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PhylogenyDesignTest(unittest.TestCase):
    def test_analysis_products_have_nonoverlapping_biological_contracts(self) -> None:
        rows = read_tsv(CONFIG / "analysis_designs.tsv")
        designs = {row["design_id"]: row for row in rows}
        self.assertEqual(len(designs), len(rows))
        self.assertEqual(
            set(designs),
            {
                "species_primary",
                "assembly_unit_diagnostic",
                "species_representative_swap",
                "legacy_manuscript_reproduction",
            },
        )

        primary = designs["species_primary"]
        self.assertEqual(primary["analysis_level"], "selected_species_or_parental_lineage")
        self.assertIn("OrthoFinder 3", primary["orthogroup_engine"])
        self.assertIn("IQ-TREE 2", primary["concatenation_engine"])
        self.assertIn("ASTRAL-Pro", primary["coalescent_engine"])
        self.assertIn("MCMCTree", primary["dating_policy"])
        self.assertIn("CAFE 5", primary["gene_family_policy"])
        self.assertIn("Enabled only after", primary["pgls_policy"])
        self.assertIn("A. x zhejiangensis A and B", primary["pgls_policy"])
        self.assertIn("A. rufa", primary["pgls_policy"])
        self.assertIn("leave-one-lineage-out", primary["pgls_policy"])

        diagnostic = designs["assembly_unit_diagnostic"]
        self.assertEqual(diagnostic["analysis_level"], "assembly_unit")
        self.assertEqual(diagnostic["dating_policy"], "prohibited")
        self.assertEqual(diagnostic["gene_family_policy"], "prohibited")
        self.assertEqual(diagnostic["pgls_policy"], "prohibited")

        sensitivity = designs["species_representative_swap"]
        self.assertEqual(sensitivity["analysis_level"], "selected_species_or_parental_lineage")
        self.assertEqual(sensitivity["sensitivity_of"], "species_primary")
        self.assertIn("matching loss row", sensitivity["pgls_policy"])
        self.assertIn("time tree", sensitivity["pgls_policy"])

        legacy = designs["legacy_manuscript_reproduction"]
        self.assertEqual(legacy["analysis_level"], "legacy_mixed_tip_set")
        self.assertIn("OrthoFinder 2.5.5", legacy["orthogroup_engine"])
        self.assertIn("RAxML-NG", legacy["concatenation_engine"])
        self.assertIn("CAFE 4.2", legacy["gene_family_policy"])
        self.assertEqual(legacy["status"], "legacy_reproduction_only")

    def test_representation_policy_covers_every_focal_taxon(self) -> None:
        taxa = read_tsv(CONFIG / "taxa.tsv")
        focal = {
            row["taxon_id"]
            for row in taxa
            if row["taxon_id"] == "clem_scandens"
            or row["biological_species"].startswith("Actinidia ")
        }
        policies = read_tsv(CONFIG / "representation_policy.tsv")
        by_taxon = {row["taxon_id"]: row for row in policies}
        self.assertEqual(len(by_taxon), len(policies))
        self.assertEqual(set(by_taxon), focal)

        deliciosa = by_taxon["act_deliciosa"]
        self.assertIn("D", deliciosa["primary_preference"])
        self.assertIn("A-F", deliciosa["diagnostic_representation"])
        self.assertIn("never sum A-D as four species", by_taxon["act_arguta"]["cafe5_count_policy"])

        eriantha = by_taxon["act_eriantha"]
        self.assertIn("HAP1", eriantha["primary_preference"])
        self.assertIn("HAP2", eriantha["required_sensitivities"])
        self.assertIn(
            "one non-hybrid individual", eriantha["primary_representation_rule"]
        )

        rufa = by_taxon["act_rufa"]
        self.assertIn("ARU", rufa["primary_preference"])
        self.assertIn("Fuchu", rufa["primary_preference"])
        self.assertIn("independently", rufa["primary_representation_rule"])
        self.assertIn("end-to-end representative swap", rufa["required_sensitivities"])
        self.assertIn("matching loss row", rufa["required_sensitivities"])
        self.assertIn("omit-species sensitivity", rufa["required_sensitivities"])

        zhejiangensis = by_taxon["act_zhejiangensis"]
        self.assertIn("one family-count row for A", zhejiangensis["cafe5_count_policy"])
        self.assertEqual(
            zhejiangensis["selection_status"], "selected_two_parental_lineages"
        )
        taxa_by_id = {row["taxon_id"]: row for row in taxa}
        self.assertEqual(taxa_by_id["act_zhejiangensis"]["include_species_tree"], "parental_lineage_model")

    def test_workers_tools_pgls_and_cafe_cohorts_are_frozen(self) -> None:
        parameters = tomllib.loads(
            (ROOT / "config" / "analysis_parameters.toml").read_text(
                encoding="utf-8"
            )
        )
        phylogeny = parameters["phylogeny"]
        self.assertLessEqual(phylogeny["maximum_workers"], 15)
        self.assertLessEqual(
            phylogeny["maximum_workers"], parameters["maximum_project_workers"]
        )
        self.assertIn("OrthoFinder 3", phylogeny["orthogroup_engine"])
        self.assertIn("IQ-TREE 2", phylogeny["concatenation_engine"])
        self.assertIn("ASTRAL-Pro", phylogeny["coalescent_engine"])
        self.assertIn("CAFE 5", phylogeny["gene_family_engine"])
        self.assertFalse(phylogeny["assembly_unit_dating_allowed"])
        self.assertFalse(phylogeny["assembly_unit_cafe_allowed"])
        self.assertFalse(phylogeny["polyploid_raw_unit_sum_allowed"])
        self.assertTrue(phylogeny["active_fossil_requires_exact_bracketing_taxa"])
        self.assertTrue(phylogeny["pgls"]["enabled"])
        self.assertEqual(phylogeny["pgls"]["status"], "enabled")
        self.assertIn("dated primary lineage tree PASS", phylogeny["pgls"]["execution_gate"])
        self.assertIn(
            "shared/non-shared lineage loss table PASS",
            phylogeny["pgls"]["execution_gate"],
        )
        self.assertEqual(phylogeny["pgls"]["primary_predictor"], "log2_ploidy")
        self.assertTrue(phylogeny["pgls"]["assembly_units_allowed"])
        self.assertFalse(phylogeny["pgls"]["shared_losses_allowed"])
        self.assertTrue(phylogeny["pgls"]["ordinary_gaussian_pgls_publication_allowed"])
        self.assertEqual(phylogeny["pgls"]["denominator_aware_model_status"], "additional_analysis")
        self.assertIn(
            "lineage-specific/non-shared positive loss count per selected species or parental lineage",
            phylogeny["pgls"]["response_numerator"],
        )
        self.assertIn(
            "excluding shared positives, uncertain calls, and non-callable genes",
            phylogeny["pgls"]["response_denominator"],
        )
        self.assertEqual(
            set(phylogeny["pgls"]["required_sensitivities"]),
            {
                "exclude Actinidia rufa",
                "swap the accepted Actinidia rufa assembly and rebuild its matched tree/loss row",
                "leave one biological species out and prune the topology for every fit",
            },
        )

        cohorts = {
            row["cohort_id"]: row for row in read_tsv(ROOT / "config" / "cohorts.tsv")
        }
        self.assertEqual(cohorts["species_tree"]["level"], "selected_species_or_parental_lineage")
        self.assertEqual(cohorts["assembly_unit_tree"]["level"], "assembly_unit")
        self.assertEqual(cohorts["cafe5_species"]["level"], "selected_species_or_parental_lineage")
        self.assertEqual(
            cohorts["pgls_species"]["status"],
            "complete_exploratory_publication_blocked_denominator_aware_model",
        )
        self.assertIn("separate A. x zhejiangensis A and B rows", cohorts["pgls_species"]["member_rule"])
        self.assertIn("denominator excludes shared positives", cohorts["pgls_species"]["member_rule"])
        self.assertIn("log2 ploidy", cohorts["pgls_species"]["purpose"])

    def test_pgls_documentation_keeps_species_and_denominator_contracts_visible(self) -> None:
        text = (ROOT / "docs" / "SPECIES_PGLS.md").read_text(encoding="utf-8")
        self.assertIn("completed as an exploratory analysis", text)
        self.assertIn("BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED", text)
        self.assertIn("D_i = |C_i \\ S|", text)
        self.assertIn("shared positives are excluded from both", text)
        self.assertIn("HAP1/HAP2", text)
        self.assertIn("*A. deliciosa* A--F", text)
        self.assertIn("*A. x zhejiangensis* A and B", text)
        self.assertIn("Swap the *A. rufa* assembly", text)
        self.assertIn("Topology-pruned leave-one-species-out", text)
        self.assertNotIn("disabled by author decision", text)

    def test_primary_lineage_sequence_manifest_has_frozen_17_tip_contract(self) -> None:
        rows = read_tsv(CONFIG / "primary_lineage_sequence_pairs.tsv")
        self.assertEqual(len(rows), 17)
        terminal_ids = {row["terminal_id"] for row in rows}
        self.assertEqual(len(terminal_ids), 17)
        self.assertIn("act_zhejiangensis_A", terminal_ids)
        self.assertIn("act_zhejiangensis_B", terminal_ids)
        self.assertIn("act_deliciosa_ADM_D", terminal_ids)
        self.assertIn("act_eriantha_HAP1", terminal_ids)
        self.assertEqual(
            {row["terminal_id"] for row in rows if row["role"] == "primary_minimal_outgroup"},
            {"coffea_arabica_E", "rhodo_simsii", "vitis_vinifera"},
        )
        self.assertTrue(all(not Path(row["protein_path"]).is_absolute() for row in rows))
        self.assertTrue(all(not Path(row["cds_path"]).is_absolute() for row in rows))

    def test_legacy_isoform_policy_reextracts_multi_isoform_tree_inputs(self) -> None:
        rows = read_tsv(CONFIG / "legacy_primary_isoform_policy.tsv")
        self.assertEqual(len(rows), 11)
        for row in rows:
            self.assertFalse(Path(row["genome_path"]).is_absolute())
            self.assertFalse(Path(row["gff_path"]).is_absolute())
            genes = int(row["source_genes"])
            mrnas = int(row["source_mrna_rows"])
            if mrnas > genes:
                self.assertIn(
                    row["primary_action"],
                    {
                        "reextract",
                        "reextract_completed",
                        "reuse_legacy_longest_span_verified",
                    },
                )
                self.assertIn("longest", row["selection_policy"])
        latifolia = {row["terminal_id"]: row for row in rows}["act_latifolia"]
        self.assertEqual(latifolia["primary_action"], "reextract_completed")
        self.assertIn("374 invalid coding genes", latifolia["note"])
        arguta = {row["terminal_id"]: row for row in rows}["act_arguta_C"]
        self.assertEqual(arguta["primary_action"], "reuse_legacy_longest_span_verified")
        self.assertIn("gffread disagreed for three proteins", arguta["note"])

    def test_no_fossil_is_active_without_exact_bracketing_taxa(self) -> None:
        rows = read_tsv(CONFIG / "calibrations.tsv")
        fossils = [
            row
            for row in rows
            if row["specimen_or_material"] != "not a fossil calibration"
        ]
        self.assertGreater(len(fossils), 0)
        active = [row for row in fossils if row["status"] == "active"]
        self.assertEqual(active, [])
        for row in fossils:
            self.assertTrue(row["required_taxon_sampling"].strip())
            self.assertNotEqual(
                row["activation_gate"],
                "passed_exact_bracketing_taxa_and_asset_qc",
            )
            self.assertTrue(row["status"].startswith("disabled_"))

        for row in active:
            self.assertEqual(
                row["activation_gate"],
                "passed_exact_bracketing_taxa_and_asset_qc",
            )


if __name__ == "__main__":
    unittest.main()
