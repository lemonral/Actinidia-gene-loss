from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _load(relative: str, name: str):
    script = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class NewLossFunctionalRendererTests(unittest.TestCase):
    def test_species_specific_figure_splits_go_and_kegg_categories(self) -> None:
        module = _load(
            "scripts/figures/render_species_specific_functional_enrichment.py",
            "species_specific_figure",
        )
        self.assertEqual(
            [label for label, _ in module.CATEGORY_COLUMNS],
            [
                "GO biological process",
                "GO molecular function",
                "GO cellular component",
                "KEGG orthology",
                "KEGG pathway",
            ],
        )

    def test_visible_branch_legend_avoids_internal_method_wording(self) -> None:
        text = (
            ROOT
            / "scripts"
            / "figures"
            / "render_tree_branch_loss_evidence.py"
        ).read_text(encoding="utf-8")
        self.assertIn('label="Branch loss"', text)
        self.assertIn('label="Strict disruption-supported subset"', text)
        self.assertNotIn('label="Article-method branch loss"', text)


if __name__ == "__main__":
    unittest.main()
