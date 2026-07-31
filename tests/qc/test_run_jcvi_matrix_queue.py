"""Small contracts for the sequential JCVI matrix materializer."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "run_jcvi_matrix_queue.py"
SPEC = importlib.util.spec_from_file_location("run_jcvi_matrix_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class MatrixQueueTests(unittest.TestCase):
    def test_under_rejects_escape(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            self.assertEqual(MODULE.under(root, "a/b"), root / "a" / "b")
            with self.assertRaises(MODULE.QueueError):
                MODULE.under(root, "../escape")

    def test_prerequisite_requires_all_units(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "state.json"
            path.write_text(
                json.dumps({"status": "PASS", "completed": [{"unit": "a"}]}) + "\n",
                encoding="utf-8",
            )
            self.assertEqual(MODULE.load_pass_state(path, "test", {"a"})["status"], "PASS")
            with self.assertRaises(MODULE.QueueError):
                MODULE.load_pass_state(path, "test", {"a", "b"})


if __name__ == "__main__":
    unittest.main()
