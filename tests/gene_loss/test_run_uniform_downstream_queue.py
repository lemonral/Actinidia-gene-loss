from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/gene_loss/run_uniform_downstream_queue.py"
SPEC = importlib.util.spec_from_file_location("uniform_downstream_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_state_writer_is_atomic(tmp_path: Path):
    path = tmp_path / "state.json"
    MODULE.write_state(path, {"status": "PASS", "rows": 817581})
    assert path.read_text().endswith("\n")
    assert '"rows": 817581' in path.read_text()
    assert list(tmp_path.iterdir()) == [path]
