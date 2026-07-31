"""Tests for the exact dated-tree to species-PGLS subtree publisher."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "prepare_species_pgls_time_tree.py"
SPEC = importlib.util.spec_from_file_location("prepare_species_pgls_time_tree", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SpeciesPGLSTimeTreeTests(unittest.TestCase):
    def fixture(self, root: Path) -> tuple[Path, Path]:
        dated = root / "dated"
        dated.mkdir()
        (dated / "dated_tree.mean_ma.tre").write_text(
            "((((A_tip:1,B_tip:1):1,(C_tip:1,D_tip:1):1):1,(E_tip:2,F_tip:2):1):2,Out:5);\n",
            encoding="utf-8",
        )
        (dated / "validation.json").write_text(json.dumps({
            "status": "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC",
            "workflow": "mcmctree_secondary_two_chain_validation_and_ultrametric_publication",
            "calibration_claim": MODULE.CALIBRATION_CLAIM,
            "dated_tree": {"tip_count": 7, "root_age_ma": 5.0},
        }) + "\n", encoding="utf-8")
        for name in ("convergence.tsv", "node_ages_ma.tsv", "secondary_calibration_summary.tsv"):
            (dated / name).write_text("fixture\n", encoding="utf-8")
        files = sorted(path for path in dated.iterdir() if path.is_file())
        with (dated / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["file", "sha256"])
            writer.writerows((path.name, digest(path)) for path in files)
        mapping = root / "map.tsv"
        mapping.write_text(
            "tree_tip\tbiological_species\tinclude\trationale\n"
            "A_tip\tSpecies A\ttrue\tprimary\n"
            "B_tip\tSpecies B\ttrue\tprimary\n"
            "C_tip\tSpecies C\ttrue\tprimary\n"
            "D_tip\tSpecies D\ttrue\tprimary\n"
            "E_tip\tSpecies E\ttrue\tprimary\n"
            "F_tip\tSpecies F\ttrue\tprimary\n"
            "Out\tOutgroup\tfalse\toutgroup\n",
            encoding="utf-8",
        )
        return dated, mapping

    def invoke(self, dated: Path, mapping: Path, output: Path) -> int:
        old = os.sys.argv
        try:
            os.sys.argv = [
                str(SCRIPT), "--dated-validation-dir", str(dated),
                "--tip-map", str(mapping), "--output-dir", str(output),
            ]
            return MODULE.main()
        finally:
            os.sys.argv = old

    def test_exact_subtree_is_atomic_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dated, mapping = self.fixture(root)
            output = root / "output"
            self.assertEqual(self.invoke(dated, mapping, output), 0)
            report = json.loads((output / "species_time_tree_pass.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(
                report["biological_species"],
                ["Species A", "Species B", "Species C", "Species D", "Species E", "Species F"],
            )
            tree = Phylo.read(str(output / "species_time_tree.nwk"), "newick")
            self.assertEqual(
                {tip.name for tip in tree.get_terminals()},
                {"Species A", "Species B", "Species C", "Species D", "Species E", "Species F"},
            )
            heights = [tree.distance(tree.root, tip) for tip in tree.get_terminals()]
            self.assertAlmostEqual(max(heights) - min(heights), 0.0)
            self.assertTrue((output / "checksums.sha256.tsv").is_file())

    def test_checksum_tampering_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dated, mapping = self.fixture(root)
            (dated / "dated_tree.mean_ma.tre").write_text("(A:1,B:1);\n", encoding="utf-8")
            output = root / "output"
            self.assertEqual(self.invoke(dated, mapping, output), 2)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
