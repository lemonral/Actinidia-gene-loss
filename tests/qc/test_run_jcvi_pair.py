"""Command-contract tests for the generic JCVI pair runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "run_jcvi_pair.py"
SPEC = importlib.util.spec_from_file_location("run_jcvi_pair", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class JcviPairCommandTest(unittest.TestCase):
    def test_commands_preserve_frozen_parameters_and_worker_count(self) -> None:
        commands = MODULE.build_commands("python3", "Clematoclethra_scandens", "Actinidia_rufa_Fuchu", 4)
        ortholog = commands["ortholog"]
        for option in (
            "--dbtype=prot",
            "--align_soft=last",
            "--cpus=4",
            "--cscore=0.7",
            "--tandem_Nmax=10",
            "--dist=20",
            "--min_size=4",
            "--no_strip_names",
            "--no_dotplot",
        ):
            self.assertIn(option, ortholog)
        self.assertEqual(commands["depth"][-1], "Clematoclethra_scandens.Actinidia_rufa_Fuchu.anchors")
        self.assertIn("--minspan=30", commands["screen"])

    def test_alias_contract_rejects_reader_facing_punctuation(self) -> None:
        with self.assertRaises(MODULE.JcviRunError):
            MODULE.validate_alias("A. rufa", "sample")
        self.assertEqual(
            MODULE.validate_alias("Actinidia_rufa_Fuchu", "sample"),
            "Actinidia_rufa_Fuchu",
        )


if __name__ == "__main__":
    unittest.main()
