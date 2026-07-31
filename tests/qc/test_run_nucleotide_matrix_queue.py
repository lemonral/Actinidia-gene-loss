from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "run_nucleotide_matrix_queue.py"
SPEC = importlib.util.spec_from_file_location("run_nucleotide_matrix_queue", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_manifest(path: Path, rows: list[str]) -> None:
    path.write_text("\t".join(MODULE.HEADER) + "\n" + "\n".join(rows) + "\n")


def test_manifest_accepts_exact_relative_paths(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    manifest = tmp_path / "manifest.tsv"
    write_manifest(
        manifest,
        ["unitA\tunitA.publisher_29chromosomes_v1\tstd/a.fa\tqc/a/target.tsv\tqc/mm/a\tqc/matrix/a"],
    )
    entries = MODULE.read_manifest(manifest, root.resolve())
    assert len(entries) == 1
    assert entries[0].unit == "unitA"
    assert entries[0].matrix_root == (root / "qc/matrix/a").resolve()


@pytest.mark.parametrize(
    "row",
    [
        "unitA\tunitA.scope\t../a.fa\tqc/a/target.tsv\tqc/mm/a\tqc/matrix/a",
        "unit A\tunitA.scope\tstd/a.fa\tqc/a/target.tsv\tqc/mm/a\tqc/matrix/a",
        "unitA\t\tstd/a.fa\tqc/a/target.tsv\tqc/mm/a\tqc/matrix/a",
    ],
)
def test_manifest_rejects_unsafe_rows(tmp_path: Path, row: str) -> None:
    root = tmp_path / "data"
    root.mkdir()
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [row])
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(manifest, root.resolve())


def test_manifest_rejects_duplicate_units_and_outputs(tmp_path: Path) -> None:
    root = tmp_path / "data"
    root.mkdir()
    manifest = tmp_path / "manifest.tsv"
    write_manifest(
        manifest,
        [
            "unitA\tunitA.scope\tstd/a.fa\tqc/a/target.tsv\tqc/mm/a\tqc/matrix/a",
            "unitA\tunitA.scope\tstd/b.fa\tqc/b/target.tsv\tqc/mm/b\tqc/matrix/b",
        ],
    )
    with pytest.raises(MODULE.QueueError):
        MODULE.read_manifest(manifest, root.resolve())


def test_capacity_gate_accepts_pass(tmp_path: Path) -> None:
    state = tmp_path / "state.json"
    state.write_text('{"status":"PASS"}\n')
    assert MODULE.wait_for_capacity_gate(state, 10)["status"] == "PASS"


def test_bundle_gate_requires_exact_unit(tmp_path: Path) -> None:
    root = tmp_path / "data"
    bundle_dir = root / "mm"
    bundle_dir.mkdir(parents=True)
    (bundle_dir / "bundle_validation.json").write_text(
        '{"status":"PASS","workflow":"bidirectional_chromosome_minimap_bundle_validation","unit":"wrong"}\n'
    )
    entry = MODULE.Entry("right", "right.scope", root / "g.fa", root / "r.tsv", bundle_dir, root / "m")
    with pytest.raises(MODULE.QueueError):
        MODULE.wait_for_bundle(entry, 10)
