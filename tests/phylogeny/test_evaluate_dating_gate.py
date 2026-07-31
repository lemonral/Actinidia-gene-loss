"""Tests for the fail-closed species-tree dating activation gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "evaluate_dating_gate.py"
SPEC = importlib.util.spec_from_file_location("evaluate_dating_gate", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class DatingGateTests(unittest.TestCase):
    def test_active_status_is_explicit(self) -> None:
        self.assertTrue(MODULE.active_status("PASS"))
        self.assertTrue(MODULE.active_status("active_primary"))
        self.assertFalse(MODULE.active_status("enabled_design_but_blocked"))
        self.assertFalse(MODULE.active_status("disabled_missing_bracket"))

    def test_binding_uses_basename_and_checksum(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.tsv"
            path.write_text("a\tb\n", encoding="utf-8")
            observed = MODULE.binding(path)
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(observed["basename"], "input.tsv")
            self.assertEqual(observed["sha256"], expected)
            self.assertNotIn(str(path.parent), json.dumps(observed))

    def test_user_authorized_secondary_status_is_active(self) -> None:
        self.assertTrue(MODULE.active_status("active_user_authorized_secondary_calibration"))


if __name__ == "__main__":
    unittest.main()
