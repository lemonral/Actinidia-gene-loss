from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "run_bidirectional_minimap_queue.py"
SPEC = importlib.util.spec_from_file_location("run_bidirectional_minimap_queue", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_manifest(path: Path, rows: list[str]) -> None:
    path.write_text("\t".join(MODULE.HEADER) + "\n" + "\n".join(rows) + "\n")


def test_manifest_accepts_prerequisite_and_two_lanes(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    for name in ("h1.fa", "a.fa", "b.fa"):
        (root / name).write_text(">Chr01\nA\n")
    manifest = tmp_path / "queue.tsv"
    write_manifest(
        manifest,
        [
            "prerequisite\t\th1\th1.fa\tqc/h1\tqc/h1.pid",
            "queued\tlane1\ta\ta.fa\tqc/a\t",
            "queued\tlane2\tb\tb.fa\tqc/b\t",
        ],
    )
    entries = MODULE.read_manifest(manifest, root.resolve())
    assert [entry.unit for entry in entries] == ["h1", "a", "b"]
    assert entries[0].controller_pid == (root / "qc" / "h1.pid").resolve()
    assert entries[1].lane == "lane1"


@pytest.mark.parametrize(
    "bad_row",
    [
        "queued\t\ta\ta.fa\tqc/a\t",
        "prerequisite\tlane1\th1\th1.fa\tqc/h1\tqc/h1.pid",
        "queued\tlane1\ta\t../a.fa\tqc/a\t",
    ],
)
def test_manifest_rejects_unsafe_rows(tmp_path: Path, bad_row: str) -> None:
    root = tmp_path / "data"
    root.mkdir()
    (root / "a.fa").write_text(">Chr01\nA\n")
    (root / "h1.fa").write_text(">Chr01\nA\n")
    manifest = tmp_path / "queue.tsv"
    write_manifest(manifest, [bad_row])
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(manifest, root.resolve())


def test_manifest_rejects_duplicate_units(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    (root / "a.fa").write_text(">Chr01\nA\n")
    manifest = tmp_path / "queue.tsv"
    write_manifest(
        manifest,
        [
            "queued\tlane1\ta\ta.fa\tqc/a\t",
            "queued\tlane2\ta\ta.fa\tqc/b\t",
        ],
    )
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(manifest, root.resolve())
