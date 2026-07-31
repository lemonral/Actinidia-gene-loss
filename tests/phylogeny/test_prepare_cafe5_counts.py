"""Unit tests for topology-bound CAFE5 count preparation."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "prepare_cafe5_counts.py"
SPEC = importlib.util.spec_from_file_location("prepare_cafe5_counts", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class CafeCountTests(unittest.TestCase):
    def test_membership_count(self) -> None:
        self.assertEqual(MODULE.membership_count(""), 0)
        self.assertEqual(MODULE.membership_count("gene1"), 1)
        self.assertEqual(MODULE.membership_count("gene1, gene2, gene3"), 3)
        with self.assertRaises(MODULE.CountError):
            MODULE.membership_count("gene1, gene1")
        with self.assertRaises(MODULE.CountError):
            MODULE.membership_count("gene1,")

    def test_terminal_manifest_requires_confirmed_unique_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "terminals.tsv"
            path.write_text(
                "terminal_id\tcanonical_tree_label\tinclude_species_tree\tidentity_status\n"
                "a\tA\ttrue\tconfirmed\n"
                "b\tB\tfalse\tpending\n",
                encoding="utf-8",
            )
            mapping, selected = MODULE.read_terminals(path)
            self.assertEqual(mapping, {"A": "a"})
            self.assertEqual(len(selected), 1)


if __name__ == "__main__":
    unittest.main()
