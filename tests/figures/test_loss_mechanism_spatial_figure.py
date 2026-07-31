"""Regression test for the loss-mechanism spatial publication figure."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "figures" / "render_loss_mechanism_spatial.py"
LOSS_TYPES = (
    "no_qualifying_translated_hit",
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
    "truncation_or_partial_alignment_candidate",
    "residual_sequence_mechanism_unresolved",
)
RELATIONS = (
    "expected_interval_local",
    "same_chromosome_displacement_candidate",
    "interchromosomal_displacement_candidate",
    "genomewide_residual_sequence_unanchored",
    "unlocalized",
)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


@unittest.skipUnless(
    importlib.util.find_spec("matplotlib") is not None
    and importlib.util.find_spec("PIL") is not None,
    "publication rendering dependencies are unavailable",
)
class LossMechanismSpatialFigureTest(unittest.TestCase):
    def test_renders_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chromosome = root / "chromosome.tsv"
            write_tsv(
                chromosome,
                [
                    {
                        "loss_type_group": loss_type,
                        "chromosome_hy4a": chrom,
                        "observed_residual_rows": "2",
                        "summed_chromosome_length_bp_across_units": "1000000",
                        "residual_rows_per_100mb": str(index + 1),
                        "length_opportunity_expected_rows": "2",
                    }
                    for index, loss_type in enumerate(LOSS_TYPES)
                    for chrom in ("Chr01", "Chr02")
                ],
            )
            location = root / "location.tsv"
            write_tsv(
                location,
                [
                    {
                        "loss_type_group": loss_type,
                        "location_relation": relation,
                        "unit_gene_rows": "2",
                    }
                    for loss_type in LOSS_TYPES
                    for relation in RELATIONS
                ],
            )
            detail = root / "detail.tsv.gz"
            detail_rows = [
                {
                    "loss_type_group": loss_type,
                    "normalized_end_distance": value,
                }
                for loss_type in LOSS_TYPES[1:4]
                for value in ("0.2", "0.5", "0.8")
            ]
            with gzip.open(detail, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(detail_rows[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(detail_rows)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "status": "PASS_LOSS_MECHANISM_SPATIAL_ANALYSIS",
                        "hy4a_standardized_chromosomes": 2,
                        "positive_unit_gene_rows": 60,
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            output = root / "figure"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--chromosome-summary",
                    str(chromosome),
                    "--location-summary",
                    str(location),
                    "--residual-positions",
                    str(detail),
                    "--run-manifest",
                    str(manifest),
                    "--output-dir",
                    str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            validation = json.loads(
                (output / "loss_mechanism_spatial.validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                validation["status"],
                "PASS_LOSS_MECHANISM_SPATIAL_FIGURE",
            )
            self.assertTrue((output / "loss_mechanism_spatial.png").is_file())


if __name__ == "__main__":
    unittest.main()
