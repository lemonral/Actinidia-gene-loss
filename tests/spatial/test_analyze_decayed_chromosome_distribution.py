"""Focused tests for decayed-only chromosome-position analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "analyze_decayed_chromosome_distribution.py"
SPEC = importlib.util.spec_from_file_location("decayed_position", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def fake_model(
    rows,
    *,
    effect,
    levels,
    adjustment_terms,
    reference_level,
):
    del rows, adjustment_terms
    summary = {
        "effect": effect,
        "model": "fixture_negative_binomial",
        "rows": 1,
        "total_decayed_loci": 1,
        "total_target_annotated_genes": 1,
        "dispersion_alpha": "0.1",
        "full_log_likelihood": "-1",
        "reduced_log_likelihood": "-2",
        "likelihood_ratio_chi_square": "2",
        "degrees_of_freedom": len(levels) - 1,
        "p_value": "0.1",
    }
    contrasts = [
        {
            effect: level,
            "comparison": (
                f"versus_{reference_level}"
                if reference_level
                else "versus_adjusted_grand_mean"
            ),
            "adjusted_decayed_loci_per_1000_genes": "1",
            "adjusted_rate_95ci_lower_per_1000": "0.5",
            "adjusted_rate_95ci_upper_per_1000": "2",
            "rate_ratio": "1",
            "rate_ratio_95ci_lower": "0.5",
            "rate_ratio_95ci_upper": "2",
            "wald_z": "0",
            "p_value": "1",
            "bh_q_value": "1",
        }
        for level in levels
    ]
    return summary, contrasts


class DecayedChromosomeDistributionTests(unittest.TestCase):
    def test_five_orientation_independent_zones(self) -> None:
        self.assertEqual(MODULE.zone_for_distance(0.0), "Z1_terminal")
        self.assertEqual(MODULE.zone_for_distance(0.2), "Z2_subterminal")
        self.assertEqual(
            MODULE.zone_for_distance(0.4),
            "Z3_intermediate_outer",
        )
        self.assertEqual(
            MODULE.zone_for_distance(0.6),
            "Z4_intermediate_inner",
        )
        self.assertEqual(MODULE.zone_for_distance(0.8), "Z5_central")
        self.assertEqual(MODULE.zone_for_distance(1.0), "Z5_central")

    def test_decayed_only_numerator_and_target_gene_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            data = root / "data"
            data.mkdir()
            registry = root / "registry.tsv"
            residual = root / "positions.tsv.gz"
            output = root / "output"
            units = ("u1", "u2")
            chromosomes = ("Chr01", "Chr02")

            registry_rows = []
            for unit in units:
                gff = data / f"{unit}.gff3"
                genome = data / f"{unit}.fa"
                with genome.open("w", encoding="utf-8") as handle:
                    for chromosome in chromosomes:
                        handle.write(f">{chromosome}\n")
                        handle.write("A" * 1000 + "\n")
                with gff.open("w", encoding="utf-8") as handle:
                    handle.write("##gff-version 3\n")
                    index = 0
                    for chromosome in chromosomes:
                        for start in (25, 225, 425, 625, 825):
                            index += 1
                            handle.write(
                                f"{chromosome}\tfixture\tgene\t{start}\t"
                                f"{start + 20}\t.\t+\t.\tID={unit}g{index}\n"
                            )
                registry_rows.append(
                    {
                        "assembly_unit_id": unit,
                        "biological_species": "Actinidia fixture",
                        "haplotype_or_subgenome": unit,
                        "source_group": "legacy" if unit == "u1" else "new",
                        "target_gff": gff.name,
                        "target_genome": genome.name,
                    }
                )
            write_tsv(
                registry,
                list(registry_rows[0]),
                registry_rows,
            )
            fields = [
                "assembly_unit_id",
                "reference_gene_id",
                "primary_classification",
                "loss_type_group",
                "residual_chromosome_hy4a",
                "residual_midpoint_1based",
                "target_chromosome_length",
                "normalized_end_distance",
                "spatial_eligible",
                "location_relation",
                "position_source",
            ]
            rows = [
                {
                    "assembly_unit_id": "u1",
                    "reference_gene_id": "r1",
                    "primary_classification": "decayed",
                    "loss_type_group": "frameshift_supported",
                    "residual_chromosome_hy4a": "Chr01",
                    "residual_midpoint_1based": "50",
                    "target_chromosome_length": "1000",
                    "normalized_end_distance": str(
                        MODULE.normalized_end_distance(50, 1000)
                    ),
                    "spatial_eligible": "true",
                    "location_relation": "expected_interval_local",
                    "position_source": "fixture",
                },
                {
                    "assembly_unit_id": "u1",
                    "reference_gene_id": "r2",
                    "primary_classification": "decayed",
                    "loss_type_group": "residual_sequence_mechanism_unresolved",
                    "residual_chromosome_hy4a": "Chr02",
                    "residual_midpoint_1based": "500",
                    "target_chromosome_length": "1000",
                    "normalized_end_distance": str(
                        MODULE.normalized_end_distance(500, 1000)
                    ),
                    "spatial_eligible": "true",
                    "location_relation": "expected_interval_local",
                    "position_source": "fixture",
                },
                {
                    "assembly_unit_id": "u2",
                    "reference_gene_id": "r3",
                    "primary_classification": "decayed",
                    "loss_type_group": "residual_sequence_mechanism_unresolved",
                    "residual_chromosome_hy4a": "",
                    "residual_midpoint_1based": "",
                    "target_chromosome_length": "",
                    "normalized_end_distance": "",
                    "spatial_eligible": "false",
                    "location_relation": "unlocalized",
                    "position_source": "",
                },
                {
                    "assembly_unit_id": "u2",
                    "reference_gene_id": "r4",
                    "primary_classification": "deleted",
                    "loss_type_group": "frameshift_supported",
                    "residual_chromosome_hy4a": "Chr01",
                    "residual_midpoint_1based": "50",
                    "target_chromosome_length": "1000",
                    "normalized_end_distance": str(
                        MODULE.normalized_end_distance(50, 1000)
                    ),
                    "spatial_eligible": "true",
                    "location_relation": "expected_interval_local",
                    "position_source": "fixture",
                },
            ]
            with gzip.open(residual, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)

            args = argparse.Namespace(
                residual_positions=residual,
                gff_registry=registry,
                data_root=data,
                output_dir=output,
                expected_units=2,
                expected_chromosomes=2,
            )
            with patch.object(
                MODULE,
                "fit_negative_binomial",
                side_effect=fake_model,
            ):
                MODULE.run(args)

            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["status"],
                "PASS_DECAYED_CHROMOSOME_POSITION_ANALYSIS",
            )
            self.assertEqual(manifest["article_method_decayed_rows"], 3)
            self.assertEqual(manifest["spatially_placed_decayed_rows"], 2)
            self.assertEqual(manifest["unlocalized_decayed_rows"], 1)
            self.assertEqual(manifest["strict_pseudogenized_placed_rows"], 1)
            self.assertEqual(manifest["target_annotated_gene_rows"], 20)

            with (
                output / "pooled_chromosome_decayed_burden.tsv"
            ).open(encoding="utf-8") as handle:
                pooled = list(csv.DictReader(handle, delimiter="\t"))
            all_decayed = [
                row for row in pooled if row["analysis_group"] == "all_decayed"
            ]
            self.assertEqual(
                sum(int(row["decayed_loci"]) for row in all_decayed),
                2,
            )
            self.assertEqual(
                sum(int(row["target_annotated_genes"]) for row in all_decayed),
                20,
            )


if __name__ == "__main__":
    unittest.main()
