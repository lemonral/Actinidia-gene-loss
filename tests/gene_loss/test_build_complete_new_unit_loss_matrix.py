"""Tests for complete callable-aware new-unit loss matrices."""

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
SCRIPT = ROOT / "scripts" / "gene_loss" / "build_complete_new_unit_loss_matrix.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def write_checksums(root: Path, names: list[str]) -> None:
    with (root / "checksums.tsv").open("w", encoding="utf-8") as handle:
        handle.write("file\tbytes\tsha256\n")
        for name in names:
            path = root / name
            handle.write(f"{name}\t{path.stat().st_size}\t{digest(path)}\n")


class CompleteNewUnitMatrixTests(unittest.TestCase):
    def _write_minimal_inputs(self, root: Path, reference: Path) -> tuple[Path, Path]:
        synorth = root / "synorth.tsv"
        synorth.write_text(
            "target1\tchr1\t1\t3\tg1\tChr01\t1\t3\t0\t0\t0\t+\tBest_Hit\n",
            encoding="utf-8",
        )
        candidates = root / "candidates"
        search = root / "search"
        candidates.mkdir()
        search.mkdir()
        candidate_table = candidates / "candidates.tsv"
        candidate_table.write_text("reference_gene\tother\ng2\tx\n", encoding="utf-8")
        (candidates / "run_manifest.json").write_text(
            json.dumps(
                {
                    "status": "PASS",
                    "unit": "u1",
                    "inputs": {
                        "reference_cds": binding(reference),
                        "synorth_pairs": binding(synorth),
                    },
                    "outputs": {"candidates": binding(candidate_table)},
                }
            )
            + "\n",
            encoding="utf-8",
        )
        columns = (
            "unit", "reference_gene", "callable", "callability_reason", "target_chromosome",
            "target_interval_start_1based", "target_interval_end_1based",
            "qualifying_genome_hit_count", "qualifying_local_hit_count", "best_hit_subject",
            "best_hit_percent_identity", "best_hit_alignment_length", "best_hit_evalue",
            "best_hit_bitscore", "primary_state", "positive_loss", "historical_reproduction_state",
        )
        with (search / "loss_states.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(
                {
                    "unit": "u1", "reference_gene": "g2", "callable": "true",
                    "callability_reason": "callable", "target_chromosome": "chr1",
                    "target_interval_start_1based": "10", "target_interval_end_1based": "20",
                    "qualifying_genome_hit_count": "0", "qualifying_local_hit_count": "0",
                    "best_hit_subject": "", "best_hit_percent_identity": "",
                    "best_hit_alignment_length": "", "best_hit_evalue": "",
                    "best_hit_bitscore": "", "primary_state": "positive_deleted",
                    "positive_loss": "true", "historical_reproduction_state": "deleted",
                }
            )
        (search / "run_manifest.json").write_text(
            json.dumps({"status": "PASS", "unit": "u1", "metrics": {"candidate_rows": 1}})
            + "\n",
            encoding="utf-8",
        )
        write_checksums(search, ["loss_states.tsv", "run_manifest.json"])
        manifest = root / "manifest.tsv"
        manifest.write_text(
            "unit\tsynorth_pairs\tcandidate_dir\tsearch_dir\n"
            "u1\tsynorth.tsv\tcandidates\tsearch\n",
            encoding="utf-8",
        )
        return manifest, search

    def test_registered_reference_symlink_may_target_legacy_store(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as legacy:
            root = Path(temporary)
            legacy_reference = Path(legacy) / "reference.fa"
            legacy_reference.write_text(">g1\nATG\n>g2\nATG\n", encoding="utf-8")
            registered_reference = root / "reference.fa"
            registered_reference.symlink_to(legacy_reference)
            manifest, _ = self._write_minimal_inputs(root, registered_reference)
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(manifest),
                    "--data-root", str(root), "--reference-cds", "reference.fa",
                    "--output-dir", str(root / "output"),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_anchor_positive_and_outside_scope_states_are_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference.fa"
            reference.write_text(">g1\nATG\n>g2\nATG\n>g3\nATG\n", encoding="utf-8")
            synorth = root / "synorth.tsv"
            synorth.write_text(
                "target1\tchr1\t1\t3\tg1\tChr01\t1\t3\t0\t0\t0\t+\tBest_Hit\n",
                encoding="utf-8",
            )
            candidates = root / "candidates"
            search = root / "search"
            candidates.mkdir()
            search.mkdir()
            candidate_table = candidates / "candidates.tsv"
            candidate_table.write_text(
                "reference_gene\tother\n"
                "g2\tx\n"
                "g_missing_cds\tx\n",
                encoding="utf-8",
            )
            candidate_manifest = {
                "status": "PASS",
                "unit": "u1",
                "inputs": {"reference_cds": binding(reference), "synorth_pairs": binding(synorth)},
                "outputs": {"candidates": binding(candidate_table)},
            }
            (candidates / "run_manifest.json").write_text(json.dumps(candidate_manifest) + "\n")

            columns = (
                "unit", "reference_gene", "callable", "callability_reason", "target_chromosome",
                "target_interval_start_1based", "target_interval_end_1based",
                "qualifying_genome_hit_count", "qualifying_local_hit_count", "best_hit_subject",
                "best_hit_percent_identity", "best_hit_alignment_length", "best_hit_evalue",
                "best_hit_bitscore", "primary_state", "positive_loss", "historical_reproduction_state",
            )
            with (search / "loss_states.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(
                    {
                        "unit": "u1", "reference_gene": "g2", "callable": "true",
                        "callability_reason": "callable", "target_chromosome": "chr1",
                        "target_interval_start_1based": "10", "target_interval_end_1based": "20",
                        "qualifying_genome_hit_count": "0", "qualifying_local_hit_count": "0",
                        "best_hit_subject": "", "best_hit_percent_identity": "",
                        "best_hit_alignment_length": "", "best_hit_evalue": "",
                        "best_hit_bitscore": "", "primary_state": "positive_deleted",
                        "positive_loss": "true", "historical_reproduction_state": "deleted",
                    }
                )
                writer.writerow(
                    {
                        "unit": "u1", "reference_gene": "g_missing_cds", "callable": "false",
                        "callability_reason": "missing_reference_cds", "target_chromosome": "",
                        "target_interval_start_1based": "", "target_interval_end_1based": "",
                        "qualifying_genome_hit_count": "0", "qualifying_local_hit_count": "0",
                        "best_hit_subject": "", "best_hit_percent_identity": "",
                        "best_hit_alignment_length": "", "best_hit_evalue": "",
                        "best_hit_bitscore": "", "primary_state": "uncertain_non_callable",
                        "positive_loss": "false", "historical_reproduction_state": "deleted",
                    }
                )
            search_manifest = {
                "status": "PASS", "unit": "u1", "metrics": {"candidate_rows": 2}
            }
            (search / "run_manifest.json").write_text(json.dumps(search_manifest) + "\n")
            write_checksums(search, ["loss_states.tsv", "run_manifest.json"])
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "unit\tsynorth_pairs\tcandidate_dir\tsearch_dir\n"
                "u1\tsynorth.tsv\tcandidates\tsearch\n",
                encoding="utf-8",
            )
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(manifest),
                    "--data-root", str(root), "--reference-cds", "reference.fa",
                    "--output-dir", str(output),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "complete_unit_loss_matrix.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                rows = {row["reference_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["g1"]["classification"], "retained")
            self.assertEqual(rows["g1"]["callable"], "true")
            self.assertEqual(rows["g2"]["classification"], "deleted")
            self.assertEqual(rows["g2"]["callable"], "true")
            self.assertEqual(rows["g3"]["classification"], "uncertain")
            self.assertEqual(rows["g3"]["callable"], "false")
            report = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(report["matrix_rows"], 3)
            self.assertEqual(report["units"][0]["candidate_rows_outside_callable_universe"], 1)
            self.assertEqual(report["definitions"]["absence_from_positive_list"], "never treated as retained")


if __name__ == "__main__":
    unittest.main()
