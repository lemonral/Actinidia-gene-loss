"""Regression test for residual-sequence loss-mechanism spatial analysis."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "analyze_loss_mechanism_positions.py"


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def paf(
    query: str,
    target: str,
    target_length: int,
    start0: int,
    end0: int,
    *,
    score: int = 200,
    frameshifts: int = 0,
    stops: int = 0,
) -> str:
    return "\t".join(
        [
            query,
            "100",
            "0",
            "100",
            "+",
            target,
            str(target_length),
            str(start0),
            str(end0),
            "90",
            "100",
            "60",
            f"AS:i:{score}",
            f"fs:i:{frameshifts}",
            f"st:i:{stops}",
        ]
    )


class LossMechanismSpatialTest(unittest.TestCase):
    def test_classifies_local_and_candidate_displacement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.tsv"
            write_tsv(
                metadata,
                [
                    {
                        "assembly_unit_id": "u1",
                        "biological_species": "Actinidia alpha",
                        "haplotype_or_subgenome": "A",
                        "include": "true",
                    },
                    {
                        "assembly_unit_id": "u2",
                        "biological_species": "Actinidia beta",
                        "haplotype_or_subgenome": "unphased",
                        "include": "true",
                    },
                ],
            )
            uniform = root / "uniform.tsv"
            write_tsv(
                uniform,
                [
                    {
                        "unit": "u1",
                        "target_genome": "genomes/u1.fa",
                        "candidate_dir": "candidates/u1",
                        "output_dir": "uniform/u1",
                    },
                    {
                        "unit": "u2",
                        "target_genome": "genomes/u2.fa",
                        "candidate_dir": "candidates/u2",
                        "output_dir": "uniform/u2",
                    },
                ],
            )
            maps = root / "maps.tsv"
            write_tsv(
                maps,
                [
                    {
                        "unit": unit,
                        "map_path": "",
                        "map_mode": "already_hy4a_chr_labels",
                    }
                    for unit in ("u1", "u2")
                ],
            )
            matrix = root / "matrix.tsv.gz"
            matrix_rows = [
                {
                    "reference_gene_id": "g1",
                    "assembly_unit_id": "u1",
                    "source_group": "legacy",
                    "manuscript_classification": "decayed",
                    "manuscript_positive_loss": "true",
                    "callable": "true",
                    "refined_decayed_cause": "frameshift_supported",
                    "refined_cause_evidence_level": "explicit_coding_disruption",
                },
                {
                    "reference_gene_id": "g2",
                    "assembly_unit_id": "u2",
                    "source_group": "new",
                    "manuscript_classification": "decayed",
                    "manuscript_positive_loss": "true",
                    "callable": "true",
                    "refined_decayed_cause": (
                        "genomewide_tblastx_hit_without_local_miniprot_support"
                    ),
                    "refined_cause_evidence_level": "sequence_detected_only",
                },
            ]
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(matrix_rows[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(matrix_rows)

            state_fields = [
                "reference_gene",
                "callable",
                "target_chromosome",
                "target_interval_start_1based",
                "target_interval_end_1based",
                "qualifying_local_alignment",
                "alignment_target_start_1based",
                "alignment_target_end_1based",
                "alignment_strand",
                "query_coverage",
                "exact_alignment_identity",
                "alignment_score",
                "frameshift_events",
                "inframe_stop_codons",
            ]
            write_tsv(
                root / "uniform" / "u1" / "uniform_candidate_loss_states.tsv",
                [
                    dict(
                        zip(
                            state_fields,
                            [
                                "g1",
                                "true",
                                "Chr01",
                                "100",
                                "500",
                                "true",
                                "200",
                                "400",
                                "+",
                                "1",
                                "0.9",
                                "200",
                                "1",
                                "0",
                            ],
                        )
                    )
                ],
            )
            write_tsv(
                root / "uniform" / "u2" / "uniform_candidate_loss_states.tsv",
                [
                    dict(
                        zip(
                            state_fields,
                            [
                                "g2",
                                "true",
                                "Chr01",
                                "100",
                                "500",
                                "false",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                            ],
                        )
                    )
                ],
            )
            for unit, lines in {
                "u1": [
                    paf("g1", "Chr01", 1000, 199, 400, frameshifts=1),
                    paf("other", "Chr02", 1200, 100, 200),
                ],
                "u2": [
                    paf("other", "Chr01", 1000, 100, 200),
                    paf("g2", "Chr02", 1200, 599, 800),
                ],
            }.items():
                with gzip.open(
                    root / "uniform" / unit / "raw_alignments.paf.gz",
                    "wt",
                    encoding="utf-8",
                ) as handle:
                    handle.write("\n".join(lines) + "\n")
            gff = root / "reference.gff3"
            gff.write_text(
                "Chr01\tsource\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                "Chr01\tsource\tmRNA\t1\t100\t.\t+\t.\tID=t1\n",
                encoding="utf-8",
            )
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--article-matrix",
                    str(matrix),
                    "--uniform-config",
                    str(uniform),
                    "--chromosome-map-config",
                    str(maps),
                    "--unit-metadata",
                    str(metadata),
                    "--reference-gff",
                    str(gff),
                    "--data-root",
                    str(root),
                    "--output-dir",
                    str(output),
                    "--expected-units",
                    "2",
                    "--expected-chromosomes",
                    "2",
                    "--expected-positive-rows",
                    "2",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["status"],
                "PASS_LOSS_MECHANISM_SPATIAL_ANALYSIS",
            )
            self.assertEqual(manifest["spatially_placed_rows"], 2)
            self.assertEqual(
                manifest["te_annotation_status"],
                "UNAVAILABLE_NO_TE_ANNOTATION",
            )
            with gzip.open(
                output / "loss_residual_positions.tsv.gz",
                "rt",
                encoding="utf-8",
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["location_relation"], "expected_interval_local")
            self.assertEqual(
                rows[1]["location_relation"],
                "interchromosomal_displacement_candidate",
            )
            self.assertEqual(rows[1]["residual_chromosome_hy4a"], "Chr02")


if __name__ == "__main__":
    unittest.main()
