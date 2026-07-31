from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "scripts"
    / "figures"
    / "render_scaffold_functional_detail.py"
)


class ScaffoldFunctionalDetailTests(unittest.TestCase):
    def test_scaffold_claims_and_counts_are_explicit(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("actual_event_gene_counts_shown", source)
        self.assertIn("point_area_linear_in_gene_count", source)
        self.assertIn("compact_bottom_layout", source)
        self.assertIn("full_term_key_in_bottom_margin", source)
        self.assertIn("resolved_23_unit_background", source)
        self.assertIn("scaffold_not_claimed_as_species_phylogeny", source)
        self.assertIn("specific_terms_named", source)
        self.assertIn("figure_title_omitted", source)
        self.assertNotIn("set_title", source)


if __name__ == "__main__":
    unittest.main()
