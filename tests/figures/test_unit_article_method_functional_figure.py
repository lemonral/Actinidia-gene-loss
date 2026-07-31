from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "scripts"
    / "figures"
    / "render_unit_article_method_functional_enrichment.py"
)


class UnitArticleMethodFunctionalFigureTests(unittest.TestCase):
    def test_publication_semantics_are_explicit(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("Lost genes (decayed + deleted)", source)
        self.assertIn("assembly_units_not_aggregated", source)
        self.assertIn("other_units_need_not_be_retained", source)
        self.assertIn("publication_title_omitted", source)
        self.assertNotIn("set_title", source)
        self.assertNotIn("Article threshold", source)


if __name__ == "__main__":
    unittest.main()
