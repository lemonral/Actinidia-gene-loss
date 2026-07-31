from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "run_jcvi_unit_queue.py"
SPEC = importlib.util.spec_from_file_location("run_jcvi_unit_queue", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def manifest(path: Path, rows: list[str]) -> None:
    path.write_text("\t".join(MODULE.HEADER) + "\n" + "\n".join(rows) + "\n")


def test_manifest_accepts_two_lanes(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    path = tmp_path / "queue.tsv"
    manifest(
        path,
        [
            "lane1\tu1\tU1\tUnit one\tacc1\tin/u1.faa\tin/u1.bed\tout/u1",
            "lane2\tu2\tU2\tUnit two\tacc2\tin/u2.faa\tin/u2.bed\tout/u2",
        ],
    )
    entries = MODULE.read_manifest(path, root.resolve())
    assert [entry.lane for entry in entries] == ["lane1", "lane2"]


@pytest.mark.parametrize(
    "row",
    [
        "\tu1\tU1\tUnit one\tacc1\tin/u1.faa\tin/u1.bed\tout/u1",
        "lane1\tu 1\tU1\tUnit one\tacc1\tin/u1.faa\tin/u1.bed\tout/u1",
        "lane1\tu1\t1U\tUnit one\tacc1\tin/u1.faa\tin/u1.bed\tout/u1",
        "lane1\tu1\tU1\tUnit one\tacc1\t../u1.faa\tin/u1.bed\tout/u1",
    ],
)
def test_manifest_rejects_unsafe_rows(tmp_path: Path, row: str) -> None:
    root = tmp_path / "data"
    root.mkdir()
    path = tmp_path / "queue.tsv"
    manifest(path, [row])
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(path, root.resolve())


def test_manifest_rejects_duplicate_alias(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    path = tmp_path / "queue.tsv"
    manifest(
        path,
        [
            "lane1\tu1\tU1\tUnit one\tacc1\tin/u1.faa\tin/u1.bed\tout/u1",
            "lane2\tu2\tU1\tUnit two\tacc2\tin/u2.faa\tin/u2.bed\tout/u2",
        ],
    )
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(path, root.resolve())
