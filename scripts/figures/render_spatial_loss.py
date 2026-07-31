#!/usr/bin/env python3
"""Render modern, opportunity-normalized gene-loss position results.

The two required inputs are ``equal_width_bins.tsv`` and
``end_distance_bins.tsv`` from ``analyze_loss_positions.py``.  The primary
panels aggregate each mutually exclusive bin across chromosomes while retaining
every assembly unit as a separate row.  Cell color is the positive-fragment
rate and every cell prints the exact numerator/denominator.

An optional ``loss_positions.tsv`` can add a descriptive centromere-distance
panel, but only for rows carrying independently supplied centromere intervals.
This panel describes the distribution among centromere-callable positive
fragments; it is explicitly not a GFF-gene-opportunity rate.  An optional
``legacy_nested_midpoint_intervals.tsv`` is rendered only as a labelled
manuscript-era sensitivity panel because its intervals overlap.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import format_taxon_label_from_metadata


PLOT_COLUMNS = (
    "panel",
    "panel_label",
    "analysis_mode",
    "analysis_label",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "assembly_scope",
    "display_label",
    "bin",
    "bin_label",
    "numerator_positive_loss_fragments",
    "denominator_count",
    "denominator_definition",
    "rate",
    "rate_percent",
    "sensitivity_only",
)

PANEL_ORDER = (
    "equal_width",
    "end_distance",
    "centromere_distance",
    "legacy_nested_sensitivity",
)

PANEL_TITLES = {
    "equal_width": "Equal-width chromosome bins",
    "end_distance": "Distance to nearest chromosome end",
    "centromere_distance": "Centromere distance (independent; descriptive)",
    "legacy_nested_sensitivity": "Legacy nested intervals (sensitivity only)",
}


class SpatialFigureError(ValueError):
    """Raised when modern spatial outputs fail publication checks."""


def _natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        raise SpatialFigureError(f"required spatial result is missing: {path.name}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or any(not field for field in fields):
            raise SpatialFigureError(f"{path.name}: missing or invalid TSV header")
        if len(fields) != len(set(fields)):
            raise SpatialFigureError(f"{path.name}: duplicate TSV column name")
        rows = list(reader)
    if not rows:
        raise SpatialFigureError(f"{path.name}: no data rows")
    return rows, fields


def _require_columns(path: Path, fields: Iterable[str], required: Iterable[str]) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise SpatialFigureError(
            f"{path.name}: missing required column(s): {', '.join(missing)}"
        )


def _parse_nonnegative_integer(value: str, *, context: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise SpatialFigureError(f"{context}: expected an integer, found {value!r}") from exc
    if parsed < 0:
        raise SpatialFigureError(f"{context}: expected a nonnegative integer")
    return parsed


def _parse_finite_float(value: str, *, context: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise SpatialFigureError(f"{context}: expected a number, found {value!r}") from exc
    if not math.isfinite(parsed):
        raise SpatialFigureError(f"{context}: expected a finite number")
    return parsed


def _validate_rate(
    value: str,
    numerator: int,
    denominator: int,
    *,
    context: str,
) -> float | None:
    if denominator == 0:
        if numerator != 0:
            raise SpatialFigureError(
                f"{context}: nonzero numerator {numerator} has zero denominator"
            )
        if value.strip():
            raise SpatialFigureError(
                f"{context}: rate must be blank when the denominator is zero"
            )
        return None
    expected = numerator / denominator
    observed = _parse_finite_float(value, context=f"{context}:reported_rate")
    if not math.isclose(observed, expected, rel_tol=1e-9, abs_tol=1e-12):
        raise SpatialFigureError(
            f"{context}: reported rate {observed} disagrees with {numerator}/{denominator}"
        )
    return expected


def _require_single_value(rows: Sequence[Mapping[str, str]], column: str, path: Path) -> str:
    values = {row[column].strip() for row in rows}
    if "" in values or len(values) != 1:
        raise SpatialFigureError(
            f"{path.name}: {column} must have exactly one nonempty value; found {sorted(values)}"
        )
    return next(iter(values))


def _unit_labels(
    metadata: Mapping[str, Mapping[str, str]],
    unit_order: Sequence[str],
) -> dict[str, str]:
    """Use an upright assembly suffix, adding scope only to resolve collisions."""

    labels = {
        unit: format_taxon_label_from_metadata(
            metadata[unit],
            suffix_fields=("haplotype_or_subgenome",),
            abbreviate_genus=True,
        )
        for unit in unit_order
    }
    collisions: dict[str, list[str]] = defaultdict(list)
    for unit, label in labels.items():
        collisions[label].append(unit)
    for duplicate_units in collisions.values():
        if len(duplicate_units) < 2:
            continue
        for unit in duplicate_units:
            labels[unit] = format_taxon_label_from_metadata(
                metadata[unit],
                suffix_fields=("haplotype_or_subgenome", "assembly_scope"),
                abbreviate_genus=True,
            )
    if len(set(labels.values())) != len(labels):
        raise SpatialFigureError(
            "reader-facing species, haplotype/subgenome, and scope metadata do not uniquely "
            "identify every assembly unit; supply distinct descriptive metadata"
        )
    return labels


def _append_plot_row(
    rows: list[dict[str, object]],
    *,
    panel: str,
    panel_label: str,
    analysis_mode: str,
    analysis_label: str,
    unit: str,
    metadata: Mapping[str, str],
    display_label: str,
    bin_number: int,
    bin_label: str,
    numerator: int,
    denominator: int,
    denominator_definition: str,
    sensitivity_only: bool,
) -> None:
    rate = numerator / denominator if denominator else None
    rows.append(
        {
            "panel": panel,
            "panel_label": panel_label,
            "analysis_mode": analysis_mode,
            "analysis_label": analysis_label,
            "assembly_unit_id": unit,
            "biological_species": metadata["biological_species"],
            "haplotype_or_subgenome": metadata["haplotype_or_subgenome"],
            "assembly_scope": metadata["assembly_scope"],
            "display_label": display_label,
            "bin": bin_number,
            "bin_label": bin_label,
            "numerator_positive_loss_fragments": numerator,
            "denominator_count": denominator,
            "denominator_definition": denominator_definition,
            "rate": rate,
            "rate_percent": 100 * rate if rate is not None else None,
            "sensitivity_only": sensitivity_only,
        }
    )


def prepare_spatial_plot(
    equal_width_bins: str | Path,
    end_distance_bins: str | Path,
    *,
    loss_positions: str | Path | None = None,
    legacy_nested_intervals: str | Path | None = None,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    """Validate modern outputs and return the exact rows used by the figure."""

    equal_path = Path(equal_width_bins)
    end_path = Path(end_distance_bins)
    equal_rows, equal_fields = _read_tsv(equal_path)
    end_rows, end_fields = _read_tsv(end_path)
    equal_required = (
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
    )
    end_required = (
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
    )
    _require_columns(equal_path, equal_fields, equal_required)
    _require_columns(end_path, end_fields, end_required)

    equal_mode = _require_single_value(equal_rows, "analysis_mode", equal_path)
    end_mode = _require_single_value(end_rows, "analysis_mode", end_path)
    if equal_mode != "primary_mutually_exclusive_equal_width":
        raise SpatialFigureError(
            f"{equal_path.name}: expected modern primary mutually exclusive equal-width bins"
        )
    if end_mode != "primary_mutually_exclusive_normalized_end_distance":
        raise SpatialFigureError(
            f"{end_path.name}: expected modern normalized end-distance bins"
        )
    equal_analysis_label = _require_single_value(equal_rows, "analysis_label", equal_path)
    end_analysis_label = _require_single_value(end_rows, "analysis_label", end_path)
    if equal_analysis_label != end_analysis_label:
        raise SpatialFigureError(
            f"{equal_path.name} and {end_path.name} have different analysis_label values"
        )
    analysis_label = equal_analysis_label

    metadata: dict[str, dict[str, str]] = {}
    end_by_unit_bin: dict[tuple[str, int], tuple[int, int]] = {}
    end_bin_bounds: dict[int, tuple[float, float]] = {}
    end_bins_by_unit: dict[str, set[int]] = defaultdict(set)
    for line_number, row in enumerate(end_rows, start=2):
        unit = row["assembly_unit_id"].strip()
        species = row["biological_species"].strip()
        suffix = row["haplotype_or_subgenome"].strip()
        scope = row["assembly_scope"].strip()
        if not all((unit, species, suffix, scope)):
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: unit and reader-facing metadata must be nonempty"
            )
        row_metadata = {
            "biological_species": species,
            "haplotype_or_subgenome": suffix,
            "assembly_scope": scope,
        }
        if unit in metadata and metadata[unit] != row_metadata:
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: inconsistent metadata for assembly unit {unit!r}"
            )
        metadata[unit] = row_metadata
        bin_number = _parse_nonnegative_integer(
            row["end_distance_bin"],
            context=f"{end_path.name}:{line_number}:end_distance_bin",
        )
        if bin_number < 1:
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: end_distance_bin must be >=1"
            )
        key = (unit, bin_number)
        if key in end_by_unit_bin:
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: duplicate unit/bin row for {unit!r}/{bin_number}"
            )
        numerator = _parse_nonnegative_integer(
            row["positive_loss_fragments"],
            context=f"{end_path.name}:{line_number}:positive_loss_fragments",
        )
        denominator = _parse_nonnegative_integer(
            row["gff_gene_opportunities"],
            context=f"{end_path.name}:{line_number}:gff_gene_opportunities",
        )
        _validate_rate(
            row["positive_loss_fragments_per_gff_gene"],
            numerator,
            denominator,
            context=f"{end_path.name}:{line_number}",
        )
        start = _parse_finite_float(
            row["normalized_end_distance_start_inclusive"],
            context=f"{end_path.name}:{line_number}:normalized_start",
        )
        finish = _parse_finite_float(
            row["normalized_end_distance_end_inclusive_only_for_last_bin"],
            context=f"{end_path.name}:{line_number}:normalized_end",
        )
        if not 0 <= start < finish <= 1:
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: normalized end-distance bounds must satisfy "
                "0 <= start < end <= 1"
            )
        bounds = (start, finish)
        if bin_number in end_bin_bounds and end_bin_bounds[bin_number] != bounds:
            raise SpatialFigureError(
                f"{end_path.name}:{line_number}: inconsistent normalized bounds for bin {bin_number}"
            )
        end_bin_bounds[bin_number] = bounds
        end_by_unit_bin[key] = (numerator, denominator)
        end_bins_by_unit[unit].add(bin_number)

    expected_bins = set(range(1, max(end_bin_bounds) + 1))
    if set(end_bin_bounds) != expected_bins:
        raise SpatialFigureError(f"{end_path.name}: end-distance bins are not contiguous from 1")
    for unit, observed in end_bins_by_unit.items():
        if observed != expected_bins:
            raise SpatialFigureError(
                f"{end_path.name}: assembly unit {unit!r} lacks the complete end-distance bin set"
            )

    equal_aggregate: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
    equal_bins_by_unit_chromosome: dict[tuple[str, str], set[int]] = defaultdict(set)
    equal_seen: set[tuple[str, str, int]] = set()
    equal_units: set[str] = set()
    for line_number, row in enumerate(equal_rows, start=2):
        unit = row["assembly_unit_id"].strip()
        chromosome = row["chromosome"].strip()
        if unit not in metadata:
            raise SpatialFigureError(
                f"{equal_path.name}:{line_number}: unit {unit!r} is absent from end_distance_bins.tsv"
            )
        if not chromosome:
            raise SpatialFigureError(f"{equal_path.name}:{line_number}: empty chromosome")
        if row["biological_species"].strip() != metadata[unit]["biological_species"] or row[
            "haplotype_or_subgenome"
        ].strip() != metadata[unit]["haplotype_or_subgenome"]:
            raise SpatialFigureError(
                f"{equal_path.name}:{line_number}: reader-facing metadata disagrees for {unit!r}"
            )
        bin_number = _parse_nonnegative_integer(
            row["bin"], context=f"{equal_path.name}:{line_number}:bin"
        )
        if bin_number not in expected_bins:
            raise SpatialFigureError(
                f"{equal_path.name}:{line_number}: bin {bin_number} is outside the end-distance bin set"
            )
        key = (unit, chromosome, bin_number)
        if key in equal_seen:
            raise SpatialFigureError(
                f"{equal_path.name}:{line_number}: duplicate unit/chromosome/bin row"
            )
        equal_seen.add(key)
        equal_bins_by_unit_chromosome[(unit, chromosome)].add(bin_number)
        numerator = _parse_nonnegative_integer(
            row["positive_loss_fragments"],
            context=f"{equal_path.name}:{line_number}:positive_loss_fragments",
        )
        denominator = _parse_nonnegative_integer(
            row["gff_gene_opportunities"],
            context=f"{equal_path.name}:{line_number}:gff_gene_opportunities",
        )
        _validate_rate(
            row["positive_loss_fragments_per_gff_gene"],
            numerator,
            denominator,
            context=f"{equal_path.name}:{line_number}",
        )
        equal_aggregate[(unit, bin_number)][0] += numerator
        equal_aggregate[(unit, bin_number)][1] += denominator
        equal_units.add(unit)

    if equal_units != set(metadata):
        raise SpatialFigureError(
            f"{equal_path.name} and {end_path.name} have different assembly-unit scopes"
        )
    for (unit, chromosome), observed in equal_bins_by_unit_chromosome.items():
        if observed != expected_bins:
            raise SpatialFigureError(
                f"{equal_path.name}: {unit!r}/{chromosome!r} lacks the complete equal-width bin set"
            )

    for unit in metadata:
        equal_total_numerator = sum(
            equal_aggregate[(unit, bin_number)][0] for bin_number in expected_bins
        )
        equal_total_denominator = sum(
            equal_aggregate[(unit, bin_number)][1] for bin_number in expected_bins
        )
        end_total_numerator = sum(
            end_by_unit_bin[(unit, bin_number)][0] for bin_number in expected_bins
        )
        end_total_denominator = sum(
            end_by_unit_bin[(unit, bin_number)][1] for bin_number in expected_bins
        )
        if (equal_total_numerator, equal_total_denominator) != (
            end_total_numerator,
            end_total_denominator,
        ):
            raise SpatialFigureError(
                f"primary panel reconciliation failed for {unit!r}: equal-width totals "
                f"{equal_total_numerator}/{equal_total_denominator} versus end-distance totals "
                f"{end_total_numerator}/{end_total_denominator}"
            )

    unit_order = sorted(
        metadata,
        key=lambda unit: (
            _natural_key(metadata[unit]["biological_species"]),
            _natural_key(metadata[unit]["haplotype_or_subgenome"]),
            _natural_key(metadata[unit]["assembly_scope"]),
            _natural_key(unit),
        ),
    )
    labels = _unit_labels(metadata, unit_order)
    number_of_bins = len(expected_bins)
    plot_rows: list[dict[str, object]] = []

    for unit in unit_order:
        for bin_number in sorted(expected_bins):
            numerator, denominator = equal_aggregate[(unit, bin_number)]
            edge_note = (
                "start" if bin_number == 1 else "end" if bin_number == number_of_bins else ""
            )
            bin_label = str(bin_number) + (f" ({edge_note})" if edge_note else "")
            _append_plot_row(
                plot_rows,
                panel="equal_width",
                panel_label="Mutually exclusive equal-width chromosome bins",
                analysis_mode=equal_mode,
                analysis_label=analysis_label,
                unit=unit,
                metadata=metadata[unit],
                display_label=labels[unit],
                bin_number=bin_number,
                bin_label=bin_label,
                numerator=numerator,
                denominator=denominator,
                denominator_definition=(
                    "GFF gene opportunities in mutually exclusive equal-width bins, summed "
                    "across gene-bearing chromosomes"
                ),
                sensitivity_only=False,
            )
        for bin_number in sorted(expected_bins):
            numerator, denominator = end_by_unit_bin[(unit, bin_number)]
            start, finish = end_bin_bounds[bin_number]
            _append_plot_row(
                plot_rows,
                panel="end_distance",
                panel_label="Normalized distance to nearest chromosome end",
                analysis_mode=end_mode,
                analysis_label=analysis_label,
                unit=unit,
                metadata=metadata[unit],
                display_label=labels[unit],
                bin_number=bin_number,
                bin_label=f"{start:g}–{finish:g}",
                numerator=numerator,
                denominator=denominator,
                denominator_definition=(
                    "GFF gene opportunities in mutually exclusive normalized nearest-end "
                    "distance bins, summed across gene-bearing chromosomes"
                ),
                sensitivity_only=False,
            )

    inputs = [equal_path, end_path]
    centromere_panel_status = "not_requested"
    if loss_positions is not None:
        positions_path = Path(loss_positions)
        position_rows, position_fields = _read_tsv(positions_path)
        _require_columns(
            positions_path,
            position_fields,
            (
                "analysis_label",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "reference_gene_id",
                "centromere_status",
                "centromere_distance_fraction_of_chromosome",
            ),
        )
        position_label = _require_single_value(position_rows, "analysis_label", positions_path)
        if position_label != analysis_label:
            raise SpatialFigureError(
                f"{positions_path.name}: analysis_label differs from primary bin inputs"
            )
        seen_position_keys: set[tuple[str, str]] = set()
        positions_per_unit: dict[str, int] = defaultdict(int)
        independent_values: dict[str, list[float]] = defaultdict(list)
        for line_number, row in enumerate(position_rows, start=2):
            unit = row["assembly_unit_id"].strip()
            gene = row["reference_gene_id"].strip()
            if unit not in metadata:
                raise SpatialFigureError(
                    f"{positions_path.name}:{line_number}: unit {unit!r} is outside primary scope"
                )
            if not gene:
                raise SpatialFigureError(
                    f"{positions_path.name}:{line_number}: empty reference_gene_id"
                )
            if row["biological_species"].strip() != metadata[unit]["biological_species"] or row[
                "haplotype_or_subgenome"
            ].strip() != metadata[unit]["haplotype_or_subgenome"]:
                raise SpatialFigureError(
                    f"{positions_path.name}:{line_number}: reader-facing metadata disagrees for {unit!r}"
                )
            key = (unit, gene)
            if key in seen_position_keys:
                raise SpatialFigureError(
                    f"{positions_path.name}:{line_number}: duplicate unit/reference-gene position"
                )
            seen_position_keys.add(key)
            positions_per_unit[unit] += 1
            status = row["centromere_status"].strip()
            fraction_text = row["centromere_distance_fraction_of_chromosome"].strip()
            if status == "independently_supplied_interval":
                fraction = _parse_finite_float(
                    fraction_text,
                    context=f"{positions_path.name}:{line_number}:centromere_distance_fraction",
                )
                if not 0 <= fraction <= 1:
                    raise SpatialFigureError(
                        f"{positions_path.name}:{line_number}: centromere distance fraction must "
                        "be in [0,1]"
                    )
                independent_values[unit].append(fraction)
            elif fraction_text:
                raise SpatialFigureError(
                    f"{positions_path.name}:{line_number}: a centromere distance is present without "
                    "an independently supplied interval"
                )
        for unit in unit_order:
            expected_position_count = sum(
                end_by_unit_bin[(unit, bin_number)][0] for bin_number in expected_bins
            )
            if positions_per_unit[unit] != expected_position_count:
                raise SpatialFigureError(
                    f"{positions_path.name}: {unit!r} contains {positions_per_unit[unit]} positions; "
                    f"expected {expected_position_count} from end_distance_bins.tsv"
                )
        independent_count = sum(len(values) for values in independent_values.values())
        if independent_count:
            centromere_panel_status = "included_descriptive_independent_intervals"
            for unit in unit_order:
                values = independent_values[unit]
                denominator = len(values)
                binned = {bin_number: 0 for bin_number in expected_bins}
                for value in values:
                    bin_number = min(
                        number_of_bins,
                        max(1, math.floor(value * number_of_bins) + 1),
                    )
                    binned[bin_number] += 1
                for bin_number in sorted(expected_bins):
                    start = (bin_number - 1) / number_of_bins
                    finish = bin_number / number_of_bins
                    _append_plot_row(
                        plot_rows,
                        panel="centromere_distance",
                        panel_label=(
                            "Descriptive distance among losses with independently mapped centromeres"
                        ),
                        analysis_mode="descriptive_positive_fragment_centromere_distance",
                        analysis_label=analysis_label,
                        unit=unit,
                        metadata=metadata[unit],
                        display_label=labels[unit],
                        bin_number=bin_number,
                        bin_label=f"{start:g}–{finish:g}",
                        numerator=binned[bin_number],
                        denominator=denominator,
                        denominator_definition=(
                            "positive loss fragments with independently supplied centromere "
                            "intervals; descriptive only, not a GFF-gene-opportunity denominator"
                        ),
                        sensitivity_only=False,
                    )
        else:
            centromere_panel_status = "omitted_no_independent_intervals"
        inputs.append(positions_path)

    legacy_panel_status = "not_requested"
    if legacy_nested_intervals is not None:
        legacy_path = Path(legacy_nested_intervals)
        legacy_rows, legacy_fields = _read_tsv(legacy_path)
        _require_columns(
            legacy_path,
            legacy_fields,
            (
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
            ),
        )
        legacy_mode = _require_single_value(legacy_rows, "analysis_mode", legacy_path)
        legacy_label = _require_single_value(legacy_rows, "analysis_label", legacy_path)
        if legacy_mode != "manuscript_era_nested_midpoint_reproduction_only":
            raise SpatialFigureError(
                f"{legacy_path.name}: legacy input is not explicitly labelled reproduction-only"
            )
        if legacy_label != analysis_label:
            raise SpatialFigureError(
                f"{legacy_path.name}: analysis_label differs from primary bin inputs"
            )
        legacy_aggregate: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])
        legacy_seen: set[tuple[str, str, int]] = set()
        legacy_bins_by_unit_chromosome: dict[tuple[str, str], set[int]] = defaultdict(set)
        legacy_units: set[str] = set()
        for line_number, row in enumerate(legacy_rows, start=2):
            unit = row["assembly_unit_id"].strip()
            chromosome = row["chromosome"].strip()
            if unit not in metadata:
                raise SpatialFigureError(
                    f"{legacy_path.name}:{line_number}: unit {unit!r} is outside primary scope"
                )
            if row["biological_species"].strip() != metadata[unit]["biological_species"] or row[
                "haplotype_or_subgenome"
            ].strip() != metadata[unit]["haplotype_or_subgenome"]:
                raise SpatialFigureError(
                    f"{legacy_path.name}:{line_number}: reader-facing metadata disagrees for {unit!r}"
                )
            if row["intervals_are_mutually_exclusive"].strip().lower() != "false" or row[
                "inferential_test_permitted"
            ].strip().lower() != "false":
                raise SpatialFigureError(
                    f"{legacy_path.name}:{line_number}: legacy intervals must be explicitly "
                    "overlapping and non-inferential"
                )
            bin_number = _parse_nonnegative_integer(
                row["nested_interval"],
                context=f"{legacy_path.name}:{line_number}:nested_interval",
            )
            if bin_number not in expected_bins:
                raise SpatialFigureError(
                    f"{legacy_path.name}:{line_number}: nested interval {bin_number} is outside "
                    "the primary bin count"
                )
            key = (unit, chromosome, bin_number)
            if key in legacy_seen:
                raise SpatialFigureError(
                    f"{legacy_path.name}:{line_number}: duplicate unit/chromosome/nested interval"
                )
            legacy_seen.add(key)
            numerator = _parse_nonnegative_integer(
                row["positive_loss_fragments"],
                context=f"{legacy_path.name}:{line_number}:positive_loss_fragments",
            )
            denominator = _parse_nonnegative_integer(
                row["gff_gene_opportunities"],
                context=f"{legacy_path.name}:{line_number}:gff_gene_opportunities",
            )
            _validate_rate(
                row["positive_loss_fragments_per_gff_gene"],
                numerator,
                denominator,
                context=f"{legacy_path.name}:{line_number}",
            )
            legacy_aggregate[(unit, bin_number)][0] += numerator
            legacy_aggregate[(unit, bin_number)][1] += denominator
            legacy_bins_by_unit_chromosome[(unit, chromosome)].add(bin_number)
            legacy_units.add(unit)
        if legacy_units != set(metadata):
            raise SpatialFigureError(
                f"{legacy_path.name}: assembly-unit scope differs from primary inputs"
            )
        for (unit, chromosome), observed in legacy_bins_by_unit_chromosome.items():
            if observed != expected_bins:
                raise SpatialFigureError(
                    f"{legacy_path.name}: {unit!r}/{chromosome!r} lacks the complete nested interval set"
                )
        for unit in unit_order:
            for bin_number in sorted(expected_bins):
                numerator, denominator = legacy_aggregate[(unit, bin_number)]
                _append_plot_row(
                    plot_rows,
                    panel="legacy_nested_sensitivity",
                    panel_label="Sensitivity only: overlapping nested intervals",
                    analysis_mode=legacy_mode,
                    analysis_label=analysis_label,
                    unit=unit,
                    metadata=metadata[unit],
                    display_label=labels[unit],
                    bin_number=bin_number,
                    bin_label=str(bin_number),
                    numerator=numerator,
                    denominator=denominator,
                    denominator_definition=(
                        "GFF genes fully contained in overlapping nested intervals; "
                        "categories are not mutually exclusive"
                    ),
                    sensitivity_only=True,
                )
        legacy_panel_status = "included_sensitivity_only_no_inference"
        inputs.append(legacy_path)

    primary_equal_losses = sum(
        int(row["numerator_positive_loss_fragments"])
        for row in plot_rows
        if row["panel"] == "equal_width"
    )
    primary_end_losses = sum(
        int(row["numerator_positive_loss_fragments"])
        for row in plot_rows
        if row["panel"] == "end_distance"
    )
    primary_equal_opportunities = sum(
        int(row["denominator_count"])
        for row in plot_rows
        if row["panel"] == "equal_width"
    )
    primary_end_opportunities = sum(
        int(row["denominator_count"])
        for row in plot_rows
        if row["panel"] == "end_distance"
    )
    if (primary_equal_losses, primary_equal_opportunities) != (
        primary_end_losses,
        primary_end_opportunities,
    ):
        raise SpatialFigureError("global primary panel reconciliation failed")

    validation: dict[str, object] = {
        "schema_version": "1.0",
        "status": "pass",
        "analysis_label": analysis_label,
        "assembly_unit_count": len(unit_order),
        "assembly_unit_order": unit_order,
        "number_of_bins": number_of_bins,
        "primary_positive_loss_fragment_count": primary_equal_losses,
        "primary_gff_gene_opportunity_count": primary_equal_opportunities,
        "centromere_panel_status": centromere_panel_status,
        "legacy_panel_status": legacy_panel_status,
        "checks": {
            "all_assembly_units_retained": True,
            "reader_facing_labels_unique": True,
            "equal_width_bins_mutually_exclusive": True,
            "reported_rates_recomputed_from_numerators_and_denominators": True,
            "primary_panel_totals_reconciled": True,
            "legacy_intervals_used_as_primary": False,
            "centromere_panel_requires_independent_intervals": True,
        },
    }
    return plot_rows, validation, inputs


def _panel_rows(
    plot_rows: Sequence[Mapping[str, object]], panel: str
) -> list[Mapping[str, object]]:
    return [row for row in plot_rows if row["panel"] == panel]


def build_spatial_figure(
    plot_rows: list[dict[str, object]], validation: Mapping[str, object]
):
    """Build aligned heatmaps with explicit numerator/denominator cell labels."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SpatialFigureError(
            "matplotlib and numpy are required to render publication figures; install the plots extra"
        ) from exc

    panels = [panel for panel in PANEL_ORDER if _panel_rows(plot_rows, panel)]
    if panels[:2] != ["equal_width", "end_distance"]:
        raise SpatialFigureError("the two modern primary panels must be present")
    unit_order = list(validation["assembly_unit_order"])
    number_of_bins = int(validation["number_of_bins"])
    display_labels = {
        str(row["assembly_unit_id"]): str(row["display_label"])
        for row in plot_rows
        if row["panel"] == "equal_width"
    }
    figure_height = max(5.2, 0.42 * len(unit_order) + 2.2)
    figure_width = max(12.0, 5.3 * len(panels))
    figure, axes_object = plt.subplots(
        1,
        len(panels),
        figsize=(figure_width, figure_height),
        constrained_layout=True,
        squeeze=False,
    )
    axes = list(axes_object[0])
    primary_rates = [
        float(row["rate"])
        for row in plot_rows
        if row["panel"] in {"equal_width", "end_distance"} and row["rate"] is not None
    ]
    primary_vmax = max(primary_rates, default=1.0)
    if primary_vmax <= 0:
        primary_vmax = 1.0
    primary_images = []

    for panel_index, (panel, axis) in enumerate(zip(panels, axes)):
        rows = _panel_rows(plot_rows, panel)
        by_key = {
            (str(row["assembly_unit_id"]), int(row["bin"])): row for row in rows
        }
        matrix = np.full((len(unit_order), number_of_bins), np.nan, dtype=float)
        annotations: list[list[str]] = [
            ["" for _ in range(number_of_bins)] for _ in unit_order
        ]
        for unit_index, unit in enumerate(unit_order):
            for bin_number in range(1, number_of_bins + 1):
                row = by_key[(unit, bin_number)]
                rate = row["rate"]
                if rate is not None:
                    matrix[unit_index, bin_number - 1] = float(rate)
                numerator = int(row["numerator_positive_loss_fragments"])
                denominator = int(row["denominator_count"])
                annotations[unit_index][bin_number - 1] = f"{numerator}/{denominator}"

        cmap_name = (
            "viridis"
            if panel in {"equal_width", "end_distance"}
            else "magma"
            if panel == "centromere_distance"
            else "cividis"
        )
        cmap = plt.get_cmap(cmap_name).with_extremes(bad="#EEEEEE")
        if panel in {"equal_width", "end_distance"}:
            vmax = primary_vmax
        else:
            values = matrix[np.isfinite(matrix)]
            vmax = float(values.max()) if values.size else 1.0
            if vmax <= 0:
                vmax = 1.0
        image = axis.imshow(matrix, aspect="auto", interpolation="nearest", vmin=0, vmax=vmax, cmap=cmap)
        if panel in {"equal_width", "end_distance"}:
            primary_images.append(image)
        bin_labels = [
            str(by_key[(unit_order[0], bin_number)]["bin_label"])
            for bin_number in range(1, number_of_bins + 1)
        ]
        axis.set_xticks(range(number_of_bins), bin_labels, rotation=40, ha="right")
        axis.set_yticks(range(len(unit_order)))
        if panel_index == 0:
            axis.set_yticklabels([display_labels[unit] for unit in unit_order])
            axis.set_ylabel("Assembly unit")
        else:
            axis.set_yticklabels([])
        panel_letter = chr(ord("a") + panel_index)
        panel_title = PANEL_TITLES[panel]
        axis.set_title(f"({panel_letter}) {panel_title}", loc="left", fontsize=10)
        axis.set_xlabel(
            "Position bin"
            if panel == "equal_width"
            else "Normalized distance bin (0 = end, 1 = center)"
            if panel == "end_distance"
            else "Distance fraction bin"
            if panel == "centromere_distance"
            else "Nested interval number"
        )
        axis.set_xticks(
            [value - 0.5 for value in range(1, number_of_bins)], minor=True
        )
        axis.set_yticks([value - 0.5 for value in range(1, len(unit_order))], minor=True)
        axis.grid(which="minor", color="white", linewidth=0.45)
        axis.tick_params(which="minor", bottom=False, left=False)
        for unit_index in range(len(unit_order)):
            for bin_index in range(number_of_bins):
                value = matrix[unit_index, bin_index]
                text_color = "white" if math.isfinite(value) and value / vmax > 0.53 else "black"
                axis.text(
                    bin_index,
                    unit_index,
                    annotations[unit_index][bin_index],
                    ha="center",
                    va="center",
                    fontsize=5.8,
                    color=text_color,
                )
        if panel not in {"equal_width", "end_distance"}:
            colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.025)
            colorbar.set_label(
                "Fraction of centromere-callable positive losses"
                if panel == "centromere_distance"
                else "Sensitivity rate in overlapping intervals",
                fontsize=8,
            )

    if primary_images:
        colorbar = figure.colorbar(primary_images[0], ax=axes[:2], fraction=0.025, pad=0.02)
        colorbar.set_label("Positive loss fragments / GFF gene opportunities", fontsize=8)
    figure.suptitle(
        "Gene-loss position rates by biological assembly unit; each cell is numerator/denominator",
        fontsize=12,
    )
    return figure


def publish_spatial_figure(
    *,
    equal_width_bins: str | Path,
    end_distance_bins: str | Path,
    output_dir: str | Path,
    loss_positions: str | Path | None = None,
    legacy_nested_intervals: str | Path | None = None,
    basename: str = "gene_loss_spatial_distribution",
    dpi: int = 300,
) -> FigureBundle:
    plot_rows, validation, input_paths = prepare_spatial_plot(
        equal_width_bins,
        end_distance_bins,
        loss_positions=loss_positions,
        legacy_nested_intervals=legacy_nested_intervals,
    )
    figure = build_spatial_figure(plot_rows, validation)
    caption = (
        "Spatial distribution of positive gene-loss fragments across all included assembly "
        "units. Panel (a) uses mutually exclusive equal-width chromosome bins; panel (b) uses "
        "mutually exclusive normalized distance to the nearest chromosome end, where 0 is an "
        "end and 1 is the chromosome center. For both primary panels, color is the positive-loss "
        "fragment count divided by all GFF gene opportunities in the same bin, and each cell "
        "prints the exact numerator/denominator. Latin binomials are italic, whereas haplotype "
        "and subgenome suffixes are upright. If present, the independently mapped centromere "
        "panel is descriptive among centromere-callable positive fragments and is not a "
        "GFF-gene-opportunity analysis. If present, the nested-interval sensitivity panel is "
        "labelled sensitivity only; its intervals overlap and no inferential test is permitted."
    )
    try:
        return write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename=basename,
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=input_paths,
            dpi=dpi,
        )
    finally:
        try:
            import matplotlib.pyplot as plt

            plt.close(figure)
        except ImportError:  # pragma: no cover - handled during construction
            pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equal-width-bins", required=True, type=Path)
    parser.add_argument("--end-distance-bins", required=True, type=Path)
    parser.add_argument("--loss-positions", type=Path)
    parser.add_argument("--legacy-nested-intervals", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="gene_loss_spatial_distribution")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        bundle = publish_spatial_figure(
            equal_width_bins=args.equal_width_bins,
            end_distance_bins=args.end_distance_bins,
            loss_positions=args.loss_positions,
            legacy_nested_intervals=args.legacy_nested_intervals,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (OSError, SpatialFigureError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"figure_bundle\t{bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
