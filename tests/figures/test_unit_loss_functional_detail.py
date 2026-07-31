from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "scripts"
    / "figures"
    / "render_unit_loss_functional_detail.py"
)


class UnitLossFunctionalDetailTests(unittest.TestCase):
    def test_publication_and_selection_rules_are_explicit(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("actual_gene_counts_shown_by_point_area", source)
        self.assertIn("specific_terms_named", source)
        self.assertIn("kegg_pathway_names_resolved", source)
        self.assertIn("go_redundancy_reduced", source)
        self.assertIn("jaccard_cutoff", source)
        self.assertIn("study_count_min", source)
        self.assertIn("fold_enrichment_min_exclusive", source)
        self.assertIn("assembly_units_not_aggregated", source)
        self.assertIn("point_area_linear_in_gene_count", source)
        self.assertIn("axes_transposed_for_readability", source)
        self.assertIn("separate_category_figures", source)
        self.assertIn("PASS_UNIT_LOSS_FUNCTIONAL_DETAIL_COLLECTION", source)
        self.assertIn("figure_title_omitted", source)
        self.assertNotIn("set_title", source)
        self.assertNotIn("Article threshold", source)


if __name__ == "__main__":
    unittest.main()
