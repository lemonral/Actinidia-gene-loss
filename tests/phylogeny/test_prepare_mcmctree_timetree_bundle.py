"""Unit tests for the TimeTree-calibrated MCMCTree input builder."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "prepare_mcmctree_timetree_bundle.py"
SPEC = importlib.util.spec_from_file_location("prepare_mcmctree_timetree_bundle", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class PrepareMCMCTreeBundleTests(unittest.TestCase):
    def test_fasta_requires_equal_codon_length(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "tiny.fa"
            fasta.write_text(">A\nACGTAA\n>B\nACGTAA\n", encoding="utf-8")
            self.assertEqual(MODULE.read_fasta(fasta), {"A": "ACGTAA", "B": "ACGTAA"})
            fasta.write_text(">A\nACGTAA\n>B\nACGTA\n", encoding="utf-8")
            with self.assertRaises(MODULE.PreparationError):
                MODULE.read_fasta(fasta)

    def test_partitioned_phylip_splits_codon_positions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "codon.phy"
            lengths = MODULE.write_partitioned_phylip(
                output,
                {"A": "AAACCC", "B": "GGGTTT"},
                ["A", "B"],
            )
            self.assertEqual(lengths, [2, 2, 2])
            observed = output.read_text(encoding="ascii")
            self.assertIn("A  AC", observed)
            self.assertIn("B  GT", observed)

    def test_calibrations_render_on_exact_mrcas(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tree_path = Path(temporary) / "tree.tre"
            tree_path.write_text("((A,B),C);\n", encoding="ascii")
            parsed = Phylo.read(str(tree_path), "newick")
            by_name = {tip.name: tip for tip in parsed.get_terminals()}
            ab = parsed.common_ancestor(by_name["A"], by_name["B"])
            rendered = MODULE.render_calibrated_tree(
                parsed,
                {id(ab): (0.05, 0.35), id(parsed.root): (1.0, 1.2)},
            )
            self.assertEqual(rendered, "((A,B) 'B(0.05,0.35)',C) 'B(1,1.2)';")

    def test_control_freezes_relaxed_clock_and_usedata(self) -> None:
        observed = MODULE.control_text(
            seed=17,
            usedata=2,
            outfile="out.txt",
            mcmcfile="mcmc.txt",
            burnin=100,
            sampfreq=10,
            nsample=50,
            root_age_upper=1.5,
        )
        self.assertIn("usedata = 2 in.BV", observed)
        self.assertIn("clock = 2", observed)
        self.assertIn("ndata = 3", observed)
        self.assertIn("RootAge = '<1.5'", observed)


if __name__ == "__main__":
    unittest.main()
