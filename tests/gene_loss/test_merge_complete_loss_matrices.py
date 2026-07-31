"""Tests for fail-closed complete-matrix merging."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "merge_complete_loss_matrices.py"
COLUMNS = (
    "reference_gene_id", "assembly_unit_id", "classification", "callable",
    "evidence_source", "primary_search_state",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def make_source(root: Path, name: str, unit: str, genes: list[str]) -> None:
    source = root / name
    source.mkdir()
    matrix = source / "complete_unit_loss_matrix.tsv"
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            writer.writerow(
                {
                    "reference_gene_id": gene, "assembly_unit_id": unit,
                    "classification": "retained", "callable": "true",
                    "evidence_source": "test", "primary_search_state": "test",
                }
            )
    report = {
        "status": "PASS", "assembly_unit_count": 1,
        "outputs": {"complete_unit_loss_matrix": binding(matrix)},
    }
    (source / "run_manifest.json").write_text(json.dumps(report) + "\n", encoding="utf-8")
    with (source / "checksums.tsv").open("w", encoding="utf-8") as handle:
        handle.write("file\tbytes\tsha256\n")
        for filename in ("complete_unit_loss_matrix.tsv", "run_manifest.json"):
            path = source / filename
            handle.write(f"{filename}\t{path.stat().st_size}\t{digest(path)}\n")


class MergeCompleteLossMatricesTests(unittest.TestCase):
    def test_merges_disjoint_complete_sources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_source(root, "a", "u1", ["g1", "g2"])
            make_source(root, "b", "u2", ["g1", "g2"])
            sources = root / "sources.tsv"
            sources.write_text(
                "source_id\tmatrix_dir\texpected_unit_count\na\ta\t1\nb\tb\t1\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--sources", str(sources),
                    "--data-root", str(root), "--expected-total-units", "2",
                    "--output-dir", str(root / "out"),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads((root / "out" / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["matrix_rows"], 4)
            self.assertEqual(report["assembly_units"], ["u1", "u2"])

    def test_different_reference_universe_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_source(root, "a", "u1", ["g1", "g2"])
            make_source(root, "b", "u2", ["g1", "g3"])
            sources = root / "sources.tsv"
            sources.write_text(
                "source_id\tmatrix_dir\texpected_unit_count\na\ta\t1\nb\tb\t1\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--sources", str(sources),
                    "--data-root", str(root), "--expected-total-units", "2",
                    "--output-dir", str(root / "out"),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("reference-gene universe differs", completed.stderr)


if __name__ == "__main__":
    unittest.main()
