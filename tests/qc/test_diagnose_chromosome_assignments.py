"""Contracts for the non-publishing chromosome diagnostic."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "diagnose_chromosome_assignments.py"
SPEC = importlib.util.spec_from_file_location("diagnose_chromosome_assignments", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class DiagnosticTests(unittest.TestCase):
    def test_under_rejects_escape(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            self.assertEqual(MODULE.under(root, "safe/path"), root / "safe" / "path")
            with self.assertRaises(MODULE.DiagnosticError):
                MODULE.under(root, "../escape")

    def test_roles_are_exactly_four(self) -> None:
        self.assertEqual(
            set(MODULE.ROLES),
            {"nucleotide_hy4a", "jcvi_hy4a", "nucleotide_hy4p", "jcvi_hy4p"},
        )

    def test_naming_policy_is_explicitly_similarity_based(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("global_one_to_one_maximum_nucleotide_similarity", source)
        self.assertIn("absolute support is QC only", source)


if __name__ == "__main__":
    unittest.main()
