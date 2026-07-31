from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "figures" / "render_lost_nlr_structural_classes.py"


class LostNlrStructuralClassFigureTests(unittest.TestCase):
    def test_renderer_has_publication_classification_contract(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("PASS_LOST_NLR_STRUCTURAL_CLASS_FIGURE", source)
        self.assertIn("Unique shared-loss NLR genes", source)
        self.assertIn("Non-shared lost NLR genes", source)
        self.assertIn("Annotated NLR genes", source)
        self.assertIn('"TIR-NBARC-LRR"', source)
        self.assertIn('"TIR-CC-NBARC-LRR"', source)
        self.assertIn("publication_title_omitted", source)
        self.assertIn("shared_loss_unique_reference_genes_separate", source)
        self.assertIn("complete_per_unit_nlr_repertoires_shown", source)
        self.assertIn(
            "class_loss_rates_use_resolved_reference_opportunities",
            source,
        )
        self.assertIn("total_and_nonshared_loss_rates_both_shown", source)
        self.assertIn("Total loss (%)", source)
        self.assertIn("plot_columns=list(plot_rows[0])", source)
        self.assertNotIn("axis.set_title", source)
        self.assertNotIn("Article", source)
        self.assertNotIn("Manuscript", source)


if __name__ == "__main__":
    unittest.main()
