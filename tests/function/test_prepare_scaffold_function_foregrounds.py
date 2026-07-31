from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREPARE = (
    ROOT
    / "scripts"
    / "function"
    / "prepare_scaffold_function_foregrounds.py"
)
ENRICH = (
    ROOT
    / "scripts"
    / "function"
    / "run_tree_aware_manuscript_enrichment.py"
)


class ScaffoldFunctionForegroundTests(unittest.TestCase):
    def test_resolved_23_unit_background_is_explicit(self) -> None:
        prepare_source = PREPARE.read_text(encoding="utf-8")
        enrichment_source = ENRICH.read_text(encoding="utf-8")
        self.assertIn("all_23_units_resolved_article_method", prepare_source)
        self.assertIn("ambiguous_not_called", prepare_source)
        self.assertIn("topology-only assembly-unit scaffold", prepare_source)
        self.assertIn("--resolved-background-scope", enrichment_source)
        self.assertIn("PASS_UNIT_SCAFFOLD_GO_KEGG", enrichment_source)
        self.assertIn('"scaffold"', enrichment_source)


if __name__ == "__main__":
    unittest.main()
