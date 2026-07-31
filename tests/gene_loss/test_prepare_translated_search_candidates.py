"""Unit tests for callable-aware translated-search candidate preparation."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "prepare_translated_search_candidates.py"
SPEC = importlib.util.spec_from_file_location("prepare_translated_search_candidates", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class CandidateTests(unittest.TestCase):
    def test_nearest_anchor_respects_window_and_direction(self) -> None:
        genes = ["a", "b", "c", "d", "e"]
        observed = {"a", "e"}
        self.assertEqual(MODULE.nearest_anchor(genes, observed, 2, -1, 3), ("a", 2))
        self.assertEqual(MODULE.nearest_anchor(genes, observed, 2, 1, 3), ("e", 2))
        self.assertIsNone(MODULE.nearest_anchor(genes, observed, 2, -1, 1))

    def test_anchor_scope_rejects_multiple_target_chromosomes(self) -> None:
        anchor = MODULE.TargetAnchor
        self.assertEqual(MODULE.anchor_scope([anchor("q1", 10, 20), anchor("q1", 30, 40)]), ("q1", 10, 40))
        self.assertIsNone(MODULE.anchor_scope([anchor("q1", 10, 20), anchor("q2", 30, 40)]))

    def test_read_synorth_accepts_legacy_reference_column(self) -> None:
        coordinate = MODULE.Coordinate("ref1", "Chr01", 1, 9, "+")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "legacy.synorths.tsv"
            path.write_text(
                "ref1\tChr01\t10\t30\tquery1\tQueryChr\t100\t120\n",
                encoding="utf-8",
            )
            anchors = MODULE.read_synorth(
                path,
                {"ref1": coordinate},
                reference_column_1based=1,
                target_gene_column_1based=5,
            )
            self.assertEqual(anchors["ref1"], [MODULE.TargetAnchor("QueryChr", 100, 120)])

    def test_fasta_reader_rejects_duplicate_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fa"
            path.write_text(">x\nACG\n>x\nTTT\n", encoding="utf-8")
            with self.assertRaises(MODULE.CandidateError):
                MODULE.read_fasta(path)

    def test_resolve_accepts_registered_symlink_target_outside_data_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workspace = Path(tmp)
            root = workspace / "data"
            source = workspace / "frozen" / "genome.fa"
            source.parent.mkdir()
            source.write_text(">chr1\nACGT\n", encoding="utf-8")
            link = root / "legacy_linked" / "genome.fa"
            link.parent.mkdir(parents=True)
            link.symlink_to(source)
            expected = (root.resolve() / "legacy_linked/genome.fa").absolute()
            self.assertEqual(MODULE.resolve(root.resolve(), "legacy_linked/genome.fa"), expected)

    def test_resolve_rejects_parent_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            with self.assertRaises(MODULE.CandidateError):
                MODULE.resolve(root, "../outside.fa")


if __name__ == "__main__":
    unittest.main()
