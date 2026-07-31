from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (
    ROOT
    / "scripts"
    / "function"
    / "summarize_scaffold_functional_enrichment.py"
)


class ScaffoldFunctionalSummaryTests(unittest.TestCase):
    def test_summary_keeps_node_and_terminal_scope(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("PASS_UNIT_SCAFFOLD_GO_KEGG_SUMMARY", source)
        self.assertIn("node_type", source)
        self.assertIn("loss_event_gene_count", source)
        self.assertIn("resolved across all 23 assembly units", source)


if __name__ == "__main__":
    unittest.main()
