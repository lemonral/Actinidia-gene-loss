#!/usr/bin/env python3
"""Render the strict-pseudogene observed-target-locus position analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label_from_metadata


PLOT_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "taxon_label_mathtext",
    "source_group",
    "end_to_center_bin",
    "normalized_end_distance_lower",
    "normalized_end_distance_upper",
    "pseudogenized_count",
    "observed_locus_denominator",
    "pseudogenized_rate",
)


class PositionFigureError(ValueError):
    pass


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise PositionFigureError(f"missing or empty input: {path.name}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not reader.fieldnames or not rows:
        raise PositionFigureError(f"invalid or empty table: {path.name}")
    return rows


def prepare_position_plot(
    spatial_dir: str | Path,
    unit_metadata: str | Path,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    source = Path(spatial_dir)
    bin_path = source / "position_bin_summary.tsv"
    coefficients_path = source / "model_coefficients.tsv"
    positions_path = source / "pseudogenized_target_positions.tsv"
    manifest_path = source / "run_manifest.json"
    metadata_path = Path(unit_metadata)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PositionFigureError("position run manifest is missing or invalid") from exc
    if (
        manifest.get("status") != "PASS"
        or manifest.get("workflow") != "uniform_pseudogenized_target_position_analysis"
    ):
        raise PositionFigureError("position run is not the strict observed-locus PASS workflow")

    metadata_rows = read_tsv(metadata_path)
    included = [row for row in metadata_rows if row.get("include") == "true"]
    metadata = {row["assembly_unit_id"]: row for row in included}
    if len(included) != 23 or len(metadata) != 23:
        raise PositionFigureError("unit metadata must contain exactly 23 unique included units")

    bin_rows = read_tsv(bin_path)
    units = []
    seen_units = set()
    plot_rows: list[dict[str, object]] = []
    bin_counts: dict[str, set[int]] = {}
    total_pseudogenized = 0
    for line_number, row in enumerate(bin_rows, 2):
        unit = row.get("assembly_unit_id", "")
        if unit not in metadata:
            raise PositionFigureError(f"position_bin_summary.tsv:{line_number}: unknown unit")
        if unit not in seen_units:
            units.append(unit)
            seen_units.add(unit)
        try:
            bin_number = int(row["end_to_center_bin"])
            lower = float(row["normalized_end_distance_lower"])
            upper = float(row["normalized_end_distance_upper"])
            positive = int(row["primary_pseudogenized"])
            denominator = int(row["primary_observed_locus_denominator"])
            reported_rate = float(row["primary_pseudogenized_rate"])
        except (KeyError, ValueError) as exc:
            raise PositionFigureError(
                f"position_bin_summary.tsv:{line_number}: invalid primary values"
            ) from exc
        if denominator <= 0 or positive < 0 or positive > denominator:
            raise PositionFigureError(f"position_bin_summary.tsv:{line_number}: invalid count")
        expected_rate = positive / denominator
        if not math.isclose(reported_rate, expected_rate, rel_tol=0, abs_tol=1e-12):
            raise PositionFigureError(f"position_bin_summary.tsv:{line_number}: rate mismatch")
        if not (1 <= bin_number <= 10 and math.isclose(lower, (bin_number - 1) / 10)):
            raise PositionFigureError(f"position_bin_summary.tsv:{line_number}: invalid bin")
        bin_counts.setdefault(unit, set()).add(bin_number)
        total_pseudogenized += positive
        meta = metadata[unit]
        plot_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": meta["biological_species"],
                "haplotype_or_subgenome": meta["haplotype_or_subgenome"],
                "taxon_label_mathtext": format_downstream_taxon_label_from_metadata(
                    meta,
                    suffix_fields=("haplotype_or_subgenome",),
                    abbreviate_genus=True,
                ),
                "source_group": row["source_group"],
                "end_to_center_bin": bin_number,
                "normalized_end_distance_lower": lower,
                "normalized_end_distance_upper": upper,
                "pseudogenized_count": positive,
                "observed_locus_denominator": denominator,
                "pseudogenized_rate": reported_rate,
            }
        )
    if units != list(metadata) or any(value != set(range(1, 11)) for value in bin_counts.values()):
        raise PositionFigureError("bin table does not close to the exact 23-unit by 10-bin grid")
    if total_pseudogenized != 20046:
        raise PositionFigureError("strict pseudogenized total is not 20,046")

    position_rows = read_tsv(positions_path)
    if len(position_rows) != 20046:
        raise PositionFigureError("target-position table does not contain exactly 20,046 rows")
    for line_number, row in enumerate(position_rows, 2):
        try:
            start = int(row["position_start_1based"])
            end = int(row["position_end_1based"])
            frameshifts = int(row["frameshift_events"])
            stops = int(row["inframe_stop_codons"])
        except (KeyError, ValueError) as exc:
            raise PositionFigureError(
                f"pseudogenized_target_positions.tsv:{line_number}: invalid coordinate/event"
            ) from exc
        if start < 1 or end < start or frameshifts + stops <= 0:
            raise PositionFigureError(
                f"pseudogenized_target_positions.tsv:{line_number}: unsupported target locus"
            )

    coefficients = read_tsv(coefficients_path)
    selected = [
        row
        for row in coefficients
        if row.get("model") == "pooled_pseudogenized_vs_retained_unit_fe"
        and row.get("term") == "normalized_end_distance"
    ]
    if len(selected) != 1:
        raise PositionFigureError("primary pooled position coefficient is missing or duplicated")
    coefficient = selected[0]
    validation = {
        "schema_version": "1.0",
        "status": "pass",
        "analysis_role": "primary_observed_target_locus_strict_pseudogenized_vs_retained",
        "assembly_unit_count": 23,
        "bin_count_per_unit": 10,
        "primary_observed_locus_count": int(manifest["metrics"]["primary_observed_locus_rows"]),
        "strict_pseudogenized_target_locus_count": total_pseudogenized,
        "deleted_in_primary_panel": False,
        "uncertain_in_primary_panel": False,
        "pooled_position_slope_log_odds": float(coefficient["estimate_log_odds"]),
        "pooled_position_odds_ratio": float(coefficient["odds_ratio"]),
        "pooled_position_wald_p": float(coefficient["wald_p"]),
        "centromere_panel_status": "omitted_no_independent_intervals",
        "checks": {
            "actual_target_coordinates_present": True,
            "all_positions_have_strict_disruption_events": True,
            "all_23_units_present": True,
            "rates_recomputed_from_primary_counts": True,
        },
    }
    return plot_rows, validation, [bin_path, coefficients_path, positions_path, manifest_path, metadata_path]


def build_position_figure(
    plot_rows: list[dict[str, object]],
    validation: dict[str, object],
):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise PositionFigureError("matplotlib and numpy are required") from exc

    units = list(dict.fromkeys(str(row["assembly_unit_id"]) for row in plot_rows))
    labels = [
        str(next(row["taxon_label_mathtext"] for row in plot_rows if row["assembly_unit_id"] == unit))
        for unit in units
    ]
    rates = np.asarray(
        [
            [
                float(next(
                    row["pseudogenized_rate"]
                    for row in plot_rows
                    if row["assembly_unit_id"] == unit and row["end_to_center_bin"] == bin_number
                ))
                for bin_number in range(1, 11)
            ]
            for unit in units
        ]
    )
    figure, axis = plt.subplots(figsize=(11.8, 10.5), constrained_layout=True)
    image = axis.imshow(rates * 100, aspect="auto", cmap="viridis", vmin=0)
    axis.set_xticks(range(10), [f"{(i - 1) / 10:.1f}–{i / 10:.1f}" for i in range(1, 11)], rotation=35, ha="right")
    axis.set_yticks(range(len(units)), labels)
    axis.set_xlabel("Normalized distance from chromosome end (0) to center (1)")
    axis.set_ylabel("Target assembly unit")
    axis.set_title("Strict pseudogenized loci at observed target-genome positions")
    colorbar = figure.colorbar(image, ax=axis, pad=0.015)
    colorbar.set_label("Pseudogenized / observed-locus denominator (%)")
    axis.text(
        0,
        1.035,
        (
            f"Unit-fixed-effect pooled OR = {validation['pooled_position_odds_ratio']:.3f}; "
            f"clustered Wald p = {validation['pooled_position_wald_p']:.3g}; "
            "deleted and uncertain excluded"
        ),
        transform=axis.transAxes,
        fontsize=9,
        ha="left",
        va="bottom",
    )
    return figure


def publish_position_figure(
    *,
    spatial_dir: str | Path,
    unit_metadata: str | Path,
    output_dir: str | Path,
    basename: str = "uniform_pseudogenized_target_positions",
    dpi: int = 300,
) -> FigureBundle:
    rows, validation, inputs = prepare_position_plot(spatial_dir, unit_metadata)
    figure = build_position_figure(rows, validation)
    caption = (
        "Chromosome-position distribution of strict pseudogenized calls at their observed "
        "target-genome alignment loci. Each cell is the strict pseudogenized count divided "
        "by strict pseudogenized plus exact-SynOrths retained observed loci in the same "
        "assembly unit and mutually exclusive end-to-center bin. Zero is a chromosome end "
        "and one is the chromosome center. Deleted calls are absent from the primary panel "
        "because they have no observed target-gene feature; their expected-locus analysis is "
        "retained only as sensitivity evidence. Uncertain calls are excluded. Latin binomials "
        "are italic and haplotype/subgenome suffixes are upright."
    )
    try:
        return write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename=basename,
            plot_rows=rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=inputs,
            dpi=dpi,
        )
    finally:
        try:
            import matplotlib.pyplot as plt
            plt.close(figure)
        except ImportError:  # pragma: no cover
            pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial-dir", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="uniform_pseudogenized_target_positions")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        bundle = publish_position_figure(
            spatial_dir=args.spatial_dir,
            unit_metadata=args.unit_metadata,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (OSError, ValueError, TypeError, PositionFigureError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"figure_bundle\t{bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
