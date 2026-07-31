"""Focused tests for generic shared-loss and spatial publication renderers."""

from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, REPOSITORY_ROOT / "scripts" / "figures" / filename
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


species_renderer = _load_script(
    "render_species_shared_loss", "render_species_shared_loss.py"
)
spatial_renderer = _load_script("render_spatial_loss", "render_spatial_loss.py")

try:
    import matplotlib

    matplotlib.use("Agg")
except ImportError:  # pragma: no cover - optional dependency
    matplotlib = None


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


class SpeciesSharedLossRendererTest(unittest.TestCase):
    def _fixture(self, root: Path) -> Path:
        aggregate = root / "aggregate"
        aggregate.mkdir()
        species = [
            "Actinidia deliciosa",
            "Actinidia eriantha",
            "Actinidia zhejiangensis",
        ]
        statuses = {
            "g_shared": {taxon: "positive_complete" for taxon in species},
            "g_nonshared": {
                species[0]: "positive_complete",
                species[1]: "not_positive",
                species[2]: "not_positive",
            },
            "g_partial": {
                species[0]: "positive_partial",
                species[1]: "not_positive",
                species[2]: "not_positive",
            },
            "g_uncertain": {
                species[0]: "not_positive",
                species[1]: "uncertain",
                species[2]: "not_positive",
            },
            "g_none": {taxon: "not_positive" for taxon in species},
        }
        matrix_rows = []
        prevalence_rows = []
        for gene, by_species in statuses.items():
            for taxon in species:
                status = by_species[taxon]
                matrix_rows.append(
                    {
                        "reference_gene_id": gene,
                        "biological_species": taxon,
                        "species_gene_status": status,
                        "species_positive_by_rule": str(
                            status == "positive_complete"
                        ).lower(),
                    }
                )
            counts = {
                status: sum(value == status for value in by_species.values())
                for status in (
                    "positive_complete",
                    "positive_partial",
                    "not_positive",
                    "uncertain",
                )
            }
            prevalence_rows.append(
                {
                    "reference_gene_id": gene,
                    "biological_species_count": len(species),
                    "positive_complete_species_count": counts["positive_complete"],
                    "positive_partial_species_count": counts["positive_partial"],
                    "not_positive_species_count": counts["not_positive"],
                    "uncertain_species_count": counts["uncertain"],
                    "shared_positive_complete": str(
                        counts["positive_complete"] == len(species)
                    ).lower(),
                }
            )
        _write_tsv(
            aggregate / "species_gene_matrix.tsv",
            [
                "reference_gene_id",
                "biological_species",
                "species_gene_status",
                "species_positive_by_rule",
            ],
            matrix_rows,
        )
        _write_tsv(
            aggregate / "species_prevalence.tsv",
            [
                "reference_gene_id",
                "biological_species_count",
                "positive_complete_species_count",
                "positive_partial_species_count",
                "not_positive_species_count",
                "uncertain_species_count",
                "shared_positive_complete",
            ],
            prevalence_rows,
        )
        (aggregate / "species_loss_summary.json").write_text(
            json.dumps(
                {
                    "status": "complete",
                    "biological_species_count": len(species),
                    "reference_gene_count": len(statuses),
                    "shared_positive_complete_gene_count": 1,
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return aggregate

    def test_species_categories_are_disjoint_and_use_only_biological_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            aggregate = self._fixture(Path(temporary_directory))
            rows, validation, inputs = species_renderer.prepare_species_plot(aggregate)
            self.assertEqual(len(rows), 3 * 5)
            self.assertEqual(validation["biological_species_count"], 3)
            self.assertEqual(validation["reference_gene_count"], 5)
            self.assertEqual({path.name for path in inputs}, {
                "species_gene_matrix.tsv",
                "species_prevalence.tsv",
                "species_loss_summary.json",
            })

            deliciosa = {
                row["category"]: row
                for row in rows
                if row["biological_species"] == "Actinidia deliciosa"
            }
            self.assertEqual(deliciosa["shared_positive_complete"]["gene_count"], 1)
            self.assertEqual(deliciosa["non_shared_positive_complete"]["gene_count"], 1)
            self.assertEqual(deliciosa["positive_partial"]["gene_count"], 1)
            self.assertEqual(deliciosa["uncertain"]["gene_count"], 0)
            self.assertEqual(deliciosa["confidently_not_positive"]["gene_count"], 2)
            self.assertTrue(
                all(
                    sum(
                        int(row["gene_count"])
                        for row in rows
                        if row["biological_species"] == taxon
                    )
                    == 5
                    for taxon in {str(row["biological_species"]) for row in rows}
                )
            )
            self.assertTrue(
                all("technical" not in str(row["display_label"]) for row in rows)
            )
            self.assertIn(r"\mathit{A.\ deliciosa}", deliciosa["positive_partial"]["display_label"])

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_species_renderer_writes_a_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            aggregate = self._fixture(root)
            bundle = species_renderer.publish_species_figure(
                aggregate_dir=aggregate,
                output_dir=root / "figure",
                basename="shared_loss",
                dpi=72,
            )
            self.assertEqual(
                {path.name for path in bundle.directory.iterdir()},
                {
                    "shared_loss.png",
                    "shared_loss.pdf",
                    "shared_loss.plot_data.tsv",
                    "shared_loss.caption.txt",
                    "shared_loss.validation.json",
                    "shared_loss.manifest.json",
                },
            )
            self.assertIn("Partial positive evidence", bundle.plot_data.read_text())
            self.assertEqual(json.loads(bundle.validation.read_text())["status"], "pass")


class SpatialLossRendererTest(unittest.TestCase):
    UNITS = [
        *(f"deliciosa_{suffix}" for suffix in "ABCDEF"),
        *(f"zhejiangensis_{suffix}" for suffix in "AB"),
    ]

    def _fixture(self, root: Path) -> dict[str, Path]:
        equal = root / "equal_width_bins.tsv"
        end = root / "end_distance_bins.tsv"
        positions = root / "loss_positions.tsv"
        legacy = root / "legacy_nested_midpoint_intervals.tsv"
        equal_rows: list[dict[str, object]] = []
        end_rows: list[dict[str, object]] = []
        position_rows: list[dict[str, object]] = []
        legacy_rows: list[dict[str, object]] = []
        analysis_label = "primary_nonshared_pseudogenized"

        for unit_index, unit in enumerate(self.UNITS):
            if unit.startswith("deliciosa"):
                species = "Actinidia deliciosa"
            else:
                species = "Actinidia zhejiangensis"
            suffix = unit.rsplit("_", 1)[1]
            per_chromosome_losses = [unit_index % 2, 1, 0]
            if unit == self.UNITS[-1]:
                per_chromosome_losses = [0, 0, 0]
            for chromosome in ("Chr1", "Chr2"):
                for bin_number in range(1, 4):
                    numerator = per_chromosome_losses[bin_number - 1]
                    denominator = 10
                    equal_rows.append(
                        {
                            "analysis_label": analysis_label,
                            "analysis_mode": "primary_mutually_exclusive_equal_width",
                            "assembly_unit_id": unit,
                            "biological_species": species,
                            "haplotype_or_subgenome": suffix,
                            "chromosome": chromosome,
                            "bin": bin_number,
                            "gff_gene_opportunities": denominator,
                            "positive_loss_fragments": numerator,
                            "positive_loss_fragments_per_gff_gene": numerator / denominator,
                        }
                    )
                    legacy_denominator = (4, 7, 10)[bin_number - 1]
                    legacy_rows.append(
                        {
                            "analysis_label": analysis_label,
                            "analysis_mode": "manuscript_era_nested_midpoint_reproduction_only",
                            "assembly_unit_id": unit,
                            "biological_species": species,
                            "haplotype_or_subgenome": suffix,
                            "chromosome": chromosome,
                            "nested_interval": bin_number,
                            "gff_gene_opportunities": legacy_denominator,
                            "positive_loss_fragments": numerator,
                            "positive_loss_fragments_per_gff_gene": numerator / legacy_denominator,
                            "intervals_are_mutually_exclusive": "false",
                            "inferential_test_permitted": "false",
                        }
                    )
            total_positive = 2 * sum(per_chromosome_losses)
            for bin_number in range(1, 4):
                numerator = 2 * per_chromosome_losses[bin_number - 1]
                denominator = 20
                end_rows.append(
                    {
                        "analysis_label": analysis_label,
                        "analysis_mode": "primary_mutually_exclusive_normalized_end_distance",
                        "assembly_unit_id": unit,
                        "biological_species": species,
                        "haplotype_or_subgenome": suffix,
                        "assembly_scope": "chromosome_partition",
                        "end_distance_bin": bin_number,
                        "normalized_end_distance_start_inclusive": (bin_number - 1) / 3,
                        "normalized_end_distance_end_inclusive_only_for_last_bin": bin_number / 3,
                        "gff_gene_opportunities": denominator,
                        "positive_loss_fragments": numerator,
                        "positive_loss_fragments_per_gff_gene": numerator / denominator,
                    }
                )
            fractions = [0.1, 0.42, 0.88]
            for index in range(total_positive):
                position_rows.append(
                    {
                        "analysis_label": analysis_label,
                        "assembly_unit_id": unit,
                        "biological_species": species,
                        "haplotype_or_subgenome": suffix,
                        "reference_gene_id": f"{unit}_R{index + 1}",
                        "centromere_status": "independently_supplied_interval",
                        "centromere_distance_fraction_of_chromosome": fractions[index % 3],
                    }
                )

        _write_tsv(
            equal,
            [
                "analysis_label",
                "analysis_mode",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "chromosome",
                "bin",
                "gff_gene_opportunities",
                "positive_loss_fragments",
                "positive_loss_fragments_per_gff_gene",
            ],
            equal_rows,
        )
        _write_tsv(
            end,
            [
                "analysis_label",
                "analysis_mode",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "assembly_scope",
                "end_distance_bin",
                "normalized_end_distance_start_inclusive",
                "normalized_end_distance_end_inclusive_only_for_last_bin",
                "gff_gene_opportunities",
                "positive_loss_fragments",
                "positive_loss_fragments_per_gff_gene",
            ],
            end_rows,
        )
        _write_tsv(
            positions,
            [
                "analysis_label",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "reference_gene_id",
                "centromere_status",
                "centromere_distance_fraction_of_chromosome",
            ],
            position_rows,
        )
        _write_tsv(
            legacy,
            [
                "analysis_label",
                "analysis_mode",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "chromosome",
                "nested_interval",
                "gff_gene_opportunities",
                "positive_loss_fragments",
                "positive_loss_fragments_per_gff_gene",
                "intervals_are_mutually_exclusive",
                "inferential_test_permitted",
            ],
            legacy_rows,
        )
        return {
            "equal": equal,
            "end": end,
            "positions": positions,
            "legacy": legacy,
        }

    def test_all_subgenomes_rates_centromeres_and_legacy_are_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = self._fixture(Path(temporary_directory))
            rows, validation, inputs = spatial_renderer.prepare_spatial_plot(
                fixture["equal"],
                fixture["end"],
                loss_positions=fixture["positions"],
                legacy_nested_intervals=fixture["legacy"],
            )
            self.assertEqual(validation["assembly_unit_count"], 8)
            self.assertEqual(validation["number_of_bins"], 3)
            self.assertEqual(
                validation["centromere_panel_status"],
                "included_descriptive_independent_intervals",
            )
            self.assertEqual(
                validation["legacy_panel_status"],
                "included_sensitivity_only_no_inference",
            )
            self.assertEqual(len(rows), 8 * 3 * 4)
            self.assertEqual({path.name for path in inputs}, {
                "equal_width_bins.tsv",
                "end_distance_bins.tsv",
                "loss_positions.tsv",
                "legacy_nested_midpoint_intervals.tsv",
            })

            equal_a = next(
                row
                for row in rows
                if row["panel"] == "equal_width"
                and row["assembly_unit_id"] == "deliciosa_A"
                and row["bin"] == 2
            )
            self.assertEqual(equal_a["numerator_positive_loss_fragments"], 2)
            self.assertEqual(equal_a["denominator_count"], 20)
            self.assertAlmostEqual(equal_a["rate"], 0.1)
            self.assertIn("GFF gene opportunities", equal_a["denominator_definition"])

            labels = {
                row["assembly_unit_id"]: row["display_label"]
                for row in rows
                if row["panel"] == "equal_width"
            }
            self.assertIn(r"\mathit{A.\ deliciosa}", labels["deliciosa_A"])
            self.assertIn(r"\mathrm{A}", labels["deliciosa_A"])
            self.assertIn(r"\mathrm{F}", labels["deliciosa_F"])
            self.assertIn(r"\mathit{A.\ zhejiangensis}", labels["zhejiangensis_A"])
            self.assertIn(r"\mathrm{B}", labels["zhejiangensis_B"])
            self.assertTrue(
                all(
                    row["sensitivity_only"] is True
                    for row in rows
                    if row["panel"] == "legacy_nested_sensitivity"
                )
            )
            self.assertTrue(
                all(
                    row["sensitivity_only"] is False
                    for row in rows
                    if row["panel"] in {
                        "equal_width",
                        "end_distance",
                        "centromere_distance",
                    }
                )
            )

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_spatial_renderer_writes_a_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self._fixture(root)
            bundle = spatial_renderer.publish_spatial_figure(
                equal_width_bins=fixture["equal"],
                end_distance_bins=fixture["end"],
                loss_positions=fixture["positions"],
                legacy_nested_intervals=fixture["legacy"],
                output_dir=root / "figure",
                basename="spatial_loss",
                dpi=72,
            )
            self.assertTrue(bundle.png.is_file())
            self.assertTrue(bundle.pdf.is_file())
            self.assertIn("numerator_positive_loss_fragments", bundle.plot_data.read_text())
            validation = json.loads(bundle.validation.read_text())
            self.assertTrue(validation["checks"]["primary_panel_totals_reconciled"])
            self.assertFalse(validation["checks"]["legacy_intervals_used_as_primary"])

    def test_mismatched_reported_rate_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = self._fixture(Path(temporary_directory))
            text = fixture["end"].read_text(encoding="utf-8")
            fixture["end"].write_text(text.replace("\t0.1\n", "\t0.9\n", 1), encoding="utf-8")
            with self.assertRaisesRegex(
                spatial_renderer.SpatialFigureError, "reported rate"
            ):
                spatial_renderer.prepare_spatial_plot(fixture["equal"], fixture["end"])


if __name__ == "__main__":
    unittest.main()
