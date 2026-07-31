"""Small safety tests for the one-worker relabelling queue."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "run_similarity_relabel_queue.py"
SPEC = importlib.util.spec_from_file_location("run_similarity_relabel_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class QueueTests(unittest.TestCase):
    def test_resolve_rejects_escape(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            self.assertEqual(MODULE.resolve(root, "safe"), root / "safe")
            with self.assertRaises(MODULE.QueueError):
                MODULE.resolve(root, "../escape")

    def test_missing_pid_is_not_live(self) -> None:
        self.assertFalse(MODULE.live_pid(999999999))


if __name__ == "__main__":
    unittest.main()
