"""Tests for the deleted-locus spatial input adapter."""

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
SCRIPT = ROOT / "scripts" / "spatial" / "prepare_deleted_locus_spatial_inputs.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path: Path) -> dict[str, object]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def checksums(root: Path, names: list[str]) -> None:
    with (root / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("file\tbytes\tsha256\n")
        for name in names:
            path = root / name
            handle.write(f"{name}\t{path.stat().st_size}\t{digest(path)}\n")


class PrepareDeletedLocusInputsTests(unittest.TestCase):
    def test_maps_positive_expected_locus_and_binds_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            search = root / "search"
            relabel = root / "relabel"
            search.mkdir()
            relabel.mkdir()
            label = root / "labels.tsv"
            label_columns = (
                "query_chromosome", "final_chromosome", "coordinate_reference",
                "assignment_method", "assigned_score", "reciprocal_coverage",
                "orientation_to_hy4a", "hy4p_and_jcvi_agree",
                "strict_homology_gates_pass", "confidence_flag",
            )
            with label.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=label_columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for index in range(1, 30):
                    writer.writerow(
                        {
                            "query_chromosome": f"source{index}",
                            "final_chromosome": f"Chr{index:02d}",
                            "coordinate_reference": "act_chinensis_hongyang_v4_hy4a",
                            "assignment_method": "global_one_to_one_maximum_nucleotide_similarity",
                            "assigned_score": "1", "reciprocal_coverage": "1",
                            "orientation_to_hy4a": "+", "hy4p_and_jcvi_agree": "true",
                            "strict_homology_gates_pass": "true", "confidence_flag": "HIGH",
                        }
                    )

            state_columns = (
                "unit", "reference_gene", "callable", "callability_reason",
                "target_chromosome", "target_interval_start_1based",
                "target_interval_end_1based", "qualifying_genome_hit_count",
                "qualifying_local_hit_count", "best_hit_subject",
                "best_hit_percent_identity", "best_hit_alignment_length",
                "best_hit_evalue", "best_hit_bitscore", "primary_state",
                "positive_loss", "historical_reproduction_state",
            )
            with (search / "loss_states.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=state_columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(
                    {
                        "unit": "u1", "reference_gene": "g1", "callable": "true",
                        "callability_reason": "callable", "target_chromosome": "source2",
                        "target_interval_start_1based": "100", "target_interval_end_1based": "300",
                        "qualifying_genome_hit_count": "0", "qualifying_local_hit_count": "0",
                        "best_hit_subject": "", "best_hit_percent_identity": "",
                        "best_hit_alignment_length": "", "best_hit_evalue": "",
                        "best_hit_bitscore": "", "primary_state": "positive_deleted",
                        "positive_loss": "true", "historical_reproduction_state": "deleted",
                    }
                )
                writer.writerow(
                    {
                        "unit": "u1", "reference_gene": "g2", "callable": "true",
                        "callability_reason": "callable", "target_chromosome": "source3",
                        "target_interval_start_1based": "500", "target_interval_end_1based": "700",
                        "qualifying_genome_hit_count": "1", "qualifying_local_hit_count": "1",
                        "best_hit_subject": "source3", "best_hit_percent_identity": "80",
                        "best_hit_alignment_length": "100", "best_hit_evalue": "1e-20",
                        "best_hit_bitscore": "200",
                        "primary_state": "uncertain_local_genomic_sequence_detected",
                        "positive_loss": "false", "historical_reproduction_state": "decayed",
                    }
                )
            search_manifest = {
                "status": "PASS", "unit": "u1",
                "metrics": {"candidate_rows": 2, "positive_deleted": 1},
            }
            (search / "run_manifest.json").write_text(json.dumps(search_manifest) + "\n")
            checksums(search, ["loss_states.tsv", "run_manifest.json"])

            genome = relabel / "u1.genome.fa.gz"
            gff = relabel / "u1.primary.gff3"
            genome.write_bytes(b"fake-gzip-for-binding")
            gff.write_text("##gff-version 3\n", encoding="utf-8")
            relabel_manifest = {
                "status": "PASS", "unit": "u1", "inputs": {"label_map": record(label)}
            }
            (relabel / "run_manifest.json").write_text(json.dumps(relabel_manifest) + "\n")
            checksums(relabel, ["u1.genome.fa.gz", "u1.primary.gff3", "run_manifest.json"])

            manifest = root / "manifest.tsv"
            manifest.write_text(
                "unit\tbiological_species\thaplotype_or_subgenome\tassembly_scope\t"
                "search_dir\tlabel_map\trelabel_dir\n"
                "u1\tActinidia example\tA\tchromosome_partition\tsearch\tlabels.tsv\trelabel\n",
                encoding="utf-8",
            )
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--manifest", str(manifest),
                    "--data-root", str(root), "--output-dir", str(output),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "expected_deleted_locus_coordinates.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                coordinates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(coordinates), 1)
            self.assertEqual(coordinates[0]["chromosome"], "Chr02")
            self.assertEqual(coordinates[0]["expected_locus_start_1based"], "100")
            report = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["positive_deleted_count"], 1)
            self.assertIn("not an observed remnant fragment", report["coordinate_semantics"])


if __name__ == "__main__":
    unittest.main()
