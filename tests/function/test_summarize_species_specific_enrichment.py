from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "function" / "summarize_species_specific_enrichment.py"
SPEC = importlib.util.spec_from_file_location("species_specific_summary", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SpeciesSpecificSummaryTests(unittest.TestCase):
    def test_functional_categories_are_explicit(self) -> None:
        self.assertEqual(
            MODULE.category(
                {"ontology": "GO", "go_namespace": "biological_process"}
            ),
            "GO biological process",
        )
        self.assertEqual(
            MODULE.category(
                {"ontology": "GO", "go_namespace": "molecular_function"}
            ),
            "GO molecular function",
        )
        self.assertEqual(
            MODULE.category(
                {"ontology": "GO", "go_namespace": "cellular_component"}
            ),
            "GO cellular component",
        )
        self.assertEqual(
            MODULE.category({"ontology": "KEGG_KO", "go_namespace": ""}),
            "KEGG orthology",
        )
        self.assertEqual(
            MODULE.category({"ontology": "KEGG_PATHWAY", "go_namespace": ""}),
            "KEGG pathway",
        )


if __name__ == "__main__":
    unittest.main()
