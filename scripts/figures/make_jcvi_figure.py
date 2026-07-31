#!/usr/bin/env python3
"""Render a fully bidirectional JCVI coverage publication bundle.

The coverage TSV must contain one row per target ``assembly_unit_id`` and all
four required coverage dimensions: reference genes, target genes, reference
sequence bases, and target sequence bases.  Each percentage is reconciled to
its covered/total denominator before plotting.  A separate metadata TSV is
keyed by ``assembly_unit_id`` and declares each row as ``reference`` or
``target`` in ``comparison_role``.  The union of reference and target IDs in
the coverage table must exactly equal the metadata ID set; no unused or
unlabelled assembly is accepted.

Every target remains visible even when its decision status is ``excluded``.
The output is an atomic figure bundle containing PNG, PDF, exact plot data, an
English caption, JSON validation, and a checksum manifest.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import TaxonLabelError, format_taxon_label_from_metadata


SCRIPT_VERSION = "1.0.0"
DECISION_STATUSES = ("current", "candidate", "excluded")
COMPARISON_ROLES = ("reference", "target")
METADATA_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "accession",
    "assembly_scope",
    "decision_status",
    "comparison_role",
)
COVERAGE_COLUMNS = (
    "assembly_unit_id",
    "reference_assembly_unit_id",
    "decision_status",
    "reference_gene_covered",
    "reference_gene_total",
    "reference_gene_coverage_percent",
    "target_gene_covered",
    "target_gene_total",
    "target_gene_coverage_percent",
    "reference_sequence_covered_bp",
    "reference_sequence_total_bp",
    "reference_sequence_coverage_percent",
    "target_sequence_covered_bp",
    "target_sequence_total_bp",
    "target_sequence_coverage_percent",
)
PLOT_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "accession",
    "assembly_scope",
    "decision_status",
    "display_label",
    "reference_assembly_unit_id",
    "reference_biological_species",
    "reference_haplotype_or_subgenome",
    "reference_accession",
    "reference_assembly_scope",
    "reference_decision_status",
    "reference_display_label",
    "reference_gene_covered",
    "reference_gene_total",
    "reference_gene_coverage_percent",
    "target_gene_covered",
    "target_gene_total",
    "target_gene_coverage_percent",
    "reference_sequence_covered_bp",
    "reference_sequence_total_bp",
    "reference_sequence_coverage_percent",
    "target_sequence_covered_bp",
    "target_sequence_total_bp",
    "target_sequence_coverage_percent",
)
INTEGER_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)\Z")
PERCENT_TOLERANCE = 1e-6


class JCVIPlotError(RuntimeError):
    """Raised when JCVI coverage and metadata cannot be reconciled exactly."""


def read_tsv(
    path: Path, required_columns: Sequence[str], table_name: str
) -> list[dict[str, str]]:
    """Read one nonempty TSV and validate its required columns."""

    path = Path(path).expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise JCVIPlotError(f"cannot open {table_name} {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise JCVIPlotError(f"{path}: missing TSV header")
        fields = [field.strip() for field in reader.fieldnames]
        if any(not field for field in fields) or len(fields) != len(set(fields)):
            raise JCVIPlotError(f"{path}: TSV column names must be nonempty and unique")
        missing = sorted(set(required_columns).difference(fields))
        if missing:
            raise JCVIPlotError(
                f"{path}: missing {table_name} columns: {', '.join(missing)}"
            )
        reader.fieldnames = fields
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, 2):
            if None in raw:
                raise JCVIPlotError(f"{path}:{line_number}: more values than header columns")
            row = {key: (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            rows.append(row)
    if not rows:
        raise JCVIPlotError(f"{path}: {table_name} contains no data rows")
    return rows


def index_unique(
    rows: Iterable[dict[str, str]], table_name: str
) -> dict[str, dict[str, str]]:
    """Index rows by a unique, nonempty ``assembly_unit_id``."""

    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        assembly_unit_id = row["assembly_unit_id"]
        if not assembly_unit_id:
            raise JCVIPlotError(f"{table_name}: empty assembly_unit_id")
        if assembly_unit_id in indexed:
            raise JCVIPlotError(
                f"{table_name}: duplicate assembly_unit_id {assembly_unit_id!r}"
            )
        indexed[assembly_unit_id] = row
    return indexed


def parse_integer(value: str, *, field: str, assembly_unit_id: str) -> int:
    if INTEGER_PATTERN.fullmatch(value) is None:
        raise JCVIPlotError(
            f"{assembly_unit_id}: {field} must be a nonnegative integer, found {value!r}"
        )
    return int(value)


def parse_percent(value: str, *, field: str, assembly_unit_id: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise JCVIPlotError(
            f"{assembly_unit_id}: {field} must be numeric, found {value!r}"
        ) from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 100.0:
        raise JCVIPlotError(
            f"{assembly_unit_id}: {field} must be finite and within [0, 100]"
        )
    return parsed


def reconciled_metric(
    row: Mapping[str, str],
    *,
    covered_field: str,
    total_field: str,
    percent_field: str,
    assembly_unit_id: str,
) -> tuple[int, int, float]:
    """Validate one covered/total/percentage triple and return typed values."""

    covered = parse_integer(
        row[covered_field], field=covered_field, assembly_unit_id=assembly_unit_id
    )
    total = parse_integer(row[total_field], field=total_field, assembly_unit_id=assembly_unit_id)
    reported = parse_percent(
        row[percent_field], field=percent_field, assembly_unit_id=assembly_unit_id
    )
    if total == 0:
        raise JCVIPlotError(f"{assembly_unit_id}: {total_field} must be greater than zero")
    if covered > total:
        raise JCVIPlotError(
            f"{assembly_unit_id}: {covered_field} exceeds {total_field}"
        )
    expected = covered * 100.0 / total
    if abs(reported - expected) > PERCENT_TOLERANCE:
        raise JCVIPlotError(
            f"{assembly_unit_id}: {percent_field}={reported} does not reconcile to "
            f"{covered_field}/{total_field}={covered}/{total} ({expected:.9f})"
        )
    return covered, total, reported


def _label(row: Mapping[str, str], assembly_unit_id: str) -> str:
    try:
        return format_taxon_label_from_metadata(
            row,
            suffix_fields=("haplotype_or_subgenome", "accession", "decision_status"),
            abbreviate_genus=True,
        )
    except TaxonLabelError as error:
        raise JCVIPlotError(f"{assembly_unit_id}: invalid taxon metadata: {error}") from error


def prepare_plot_rows(
    metadata_path: Path, coverage_path: Path
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Validate all metadata and coverage rows and build exact plotted rows."""

    metadata_rows = read_tsv(metadata_path, METADATA_COLUMNS, "comparison metadata")
    coverage_rows = read_tsv(coverage_path, COVERAGE_COLUMNS, "JCVI coverage table")
    metadata = index_unique(metadata_rows, "comparison metadata")
    coverage = index_unique(coverage_rows, "JCVI coverage table")

    target_metadata_rows: list[dict[str, str]] = []
    reference_metadata_rows: list[dict[str, str]] = []
    for row in metadata_rows:
        assembly_unit_id = row["assembly_unit_id"]
        role = row["comparison_role"]
        if role not in COMPARISON_ROLES:
            raise JCVIPlotError(
                f"{assembly_unit_id}: comparison_role must be reference or target, found {role!r}"
            )
        if row["decision_status"] not in DECISION_STATUSES:
            raise JCVIPlotError(
                f"{assembly_unit_id}: decision_status must be one of "
                f"{', '.join(DECISION_STATUSES)}"
            )
        if not row["accession"] or not row["assembly_scope"]:
            raise JCVIPlotError(
                f"{assembly_unit_id}: accession and assembly_scope must be nonempty"
            )
        _label(row, assembly_unit_id)
        (target_metadata_rows if role == "target" else reference_metadata_rows).append(row)
    if not target_metadata_rows or not reference_metadata_rows:
        raise JCVIPlotError("comparison metadata must contain at least one reference and one target")

    target_ids = {row["assembly_unit_id"] for row in target_metadata_rows}
    coverage_target_ids = set(coverage)
    if coverage_target_ids != target_ids:
        raise JCVIPlotError(
            "JCVI target assembly_unit_id set differs from target metadata: "
            f"missing={sorted(target_ids - coverage_target_ids)}, "
            f"extra={sorted(coverage_target_ids - target_ids)}"
        )
    coverage_reference_ids = {row["reference_assembly_unit_id"] for row in coverage_rows}
    if "" in coverage_reference_ids:
        raise JCVIPlotError("JCVI coverage table contains an empty reference_assembly_unit_id")
    reference_ids = {row["assembly_unit_id"] for row in reference_metadata_rows}
    if coverage_reference_ids != reference_ids:
        raise JCVIPlotError(
            "JCVI reference assembly_unit_id set differs from reference metadata: "
            f"missing={sorted(reference_ids - coverage_reference_ids)}, "
            f"extra={sorted(coverage_reference_ids - reference_ids)}"
        )
    used_ids = coverage_target_ids.union(coverage_reference_ids)
    if used_ids != set(metadata):
        raise JCVIPlotError(
            "metadata must contain exactly the reference and target IDs used by coverage rows"
        )

    plotted: list[dict[str, object]] = []
    for plot_order, target_row in enumerate(target_metadata_rows, 1):
        assembly_unit_id = target_row["assembly_unit_id"]
        coverage_row = coverage[assembly_unit_id]
        if coverage_row["decision_status"] != target_row["decision_status"]:
            raise JCVIPlotError(
                f"{assembly_unit_id}: coverage decision_status "
                f"{coverage_row['decision_status']!r} differs from metadata "
                f"{target_row['decision_status']!r}"
            )
        reference_id = coverage_row["reference_assembly_unit_id"]
        reference_row = metadata[reference_id]

        metrics: dict[str, int | float] = {}
        for covered_field, total_field, percent_field in (
            (
                "reference_gene_covered",
                "reference_gene_total",
                "reference_gene_coverage_percent",
            ),
            ("target_gene_covered", "target_gene_total", "target_gene_coverage_percent"),
            (
                "reference_sequence_covered_bp",
                "reference_sequence_total_bp",
                "reference_sequence_coverage_percent",
            ),
            (
                "target_sequence_covered_bp",
                "target_sequence_total_bp",
                "target_sequence_coverage_percent",
            ),
        ):
            covered, total, percent = reconciled_metric(
                coverage_row,
                covered_field=covered_field,
                total_field=total_field,
                percent_field=percent_field,
                assembly_unit_id=assembly_unit_id,
            )
            metrics[covered_field] = covered
            metrics[total_field] = total
            metrics[percent_field] = percent

        plotted.append(
            {
                "plot_order": plot_order,
                "assembly_unit_id": assembly_unit_id,
                "biological_species": target_row["biological_species"],
                "haplotype_or_subgenome": target_row["haplotype_or_subgenome"],
                "accession": target_row["accession"],
                "assembly_scope": target_row["assembly_scope"],
                "decision_status": target_row["decision_status"],
                "display_label": _label(target_row, assembly_unit_id),
                "reference_assembly_unit_id": reference_id,
                "reference_biological_species": reference_row["biological_species"],
                "reference_haplotype_or_subgenome": reference_row[
                    "haplotype_or_subgenome"
                ],
                "reference_accession": reference_row["accession"],
                "reference_assembly_scope": reference_row["assembly_scope"],
                "reference_decision_status": reference_row["decision_status"],
                "reference_display_label": _label(reference_row, reference_id),
                **metrics,
            }
        )

    status_counts = Counter(str(row["decision_status"]) for row in plotted)
    validation: dict[str, object] = {
        "schema_version": "1.0",
        "renderer": "scripts/figures/make_jcvi_figure.py",
        "renderer_version": SCRIPT_VERSION,
        "status": "pass",
        "target_assembly_unit_count": len(plotted),
        "reference_assembly_unit_count": len(reference_ids),
        "excluded_target_count": status_counts.get("excluded", 0),
        "decision_status_counts": {
            status: status_counts.get(status, 0) for status in DECISION_STATUSES
        },
        "percent_reconciliation_tolerance": PERCENT_TOLERANCE,
        "checks": {
            "metadata_ids_unique": "pass",
            "coverage_target_ids_exactly_match_target_metadata": "pass",
            "coverage_reference_ids_exactly_match_reference_metadata": "pass",
            "no_unused_metadata_ids": "pass",
            "coverage_decision_status_matches_metadata": "pass",
            "all_target_rows_plotted_including_excluded": "pass",
            "all_four_bidirectional_gene_and_sequence_metrics_present": "pass",
            "all_covered_total_percent_triples_reconciled": "pass",
            "taxon_labels_metadata_driven": "pass",
        },
    }
    return plotted, validation


def make_figure(plot_rows: Sequence[Mapping[str, object]]):
    """Create four aligned JCVI coverage panels."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise JCVIPlotError(
            "Matplotlib is required to render the figure; install the project 'plots' extra"
        ) from error

    status_colors = {
        "current": "#2A9D8F",
        "candidate": "#E9C46A",
        "excluded": "#8D99AE",
    }
    status_markers = {"current": "o", "candidate": "^", "excluded": "X"}
    metric_specs = (
        ("reference_gene_coverage_percent", "Reference genes (%)"),
        ("target_gene_coverage_percent", "Target genes (%)"),
        ("reference_sequence_coverage_percent", "Reference sequence (%)"),
        ("target_sequence_coverage_percent", "Target sequence (%)"),
    )
    count = len(plot_rows)
    height = max(4.5, 0.42 * count + 1.9)
    figure, axes = plt.subplots(
        1,
        4,
        figsize=(15.5, height),
        sharey=True,
        constrained_layout=True,
    )
    positions = list(range(count))
    labels = [str(row["display_label"]) for row in plot_rows]
    for axis, (field, title) in zip(axes, metric_specs):
        axis.hlines(positions, 0, [float(row[field]) for row in plot_rows], color="#D0D0D0")
        for status in DECISION_STATUSES:
            selected = [
                (position, float(row[field]))
                for position, row in zip(positions, plot_rows)
                if row["decision_status"] == status
            ]
            if not selected:
                continue
            axis.scatter(
                [value for _, value in selected],
                [position for position, _ in selected],
                s=42,
                marker=status_markers[status],
                color=status_colors[status],
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
            )
        for position, row in zip(positions, plot_rows):
            value = float(row[field])
            axis.text(min(value + 1.0, 98.7), position, f"{value:.1f}", va="center", fontsize=7.3)
        axis.set_xlim(0, 105)
        axis.set_ylim(count - 0.5, -0.5)
        axis.set_title(title)
        axis.set_xlabel("Coverage (%)")
        axis.set_xticks((0, 25, 50, 75, 100))
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.6, alpha=0.8)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[0].set_yticks(positions, labels=labels, fontsize=8.2)
    for axis in axes[1:]:
        axis.tick_params(axis="y", labelleft=False)

    reference_labels = {str(row["reference_display_label"]) for row in plot_rows}
    if len(reference_labels) == 1:
        title = f"Bidirectional JCVI coverage against {next(iter(reference_labels))}"
    else:
        title = "Bidirectional JCVI coverage"
    figure.suptitle(title, fontsize=13)
    return figure


def caption_for(plot_rows: Sequence[Mapping[str, object]]) -> str:
    """Return an English, data-derived JCVI caption."""

    references = sorted(
        {
            f"{row['reference_biological_species']} "
            f"({row['reference_accession']}; {row['reference_assembly_scope']})"
            for row in plot_rows
        }
    )
    return (
        "Bidirectional JCVI collinearity coverage for "
        f"{len(plot_rows)} target assembly units, including targets later excluded "
        "from gene-loss analysis. The four panels separately report the percentage "
        "of eligible reference genes, target genes, reference sequence bases, and "
        "target sequence bases covered by the validated comparison. No one-directional "
        "percentage is substituted for these four measures. Latin binomials are "
        "italic; haplotype or subgenome, accession, and decision status are upright. "
        "Marker colour and shape denote current, candidate, or excluded status. "
        "Every displayed percentage is reconciled to its covered and total denominator "
        "in the plot-data table. Reference assembly units: "
        + "; ".join(references)
        + "."
    )


def render_bundle(
    *,
    metadata_path: Path,
    coverage_path: Path,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> FigureBundle:
    """Prepare, render, and atomically publish one JCVI figure bundle."""

    plot_rows, validation = prepare_plot_rows(metadata_path, coverage_path)
    figure = make_figure(plot_rows)
    try:
        return write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename=basename,
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption_for(plot_rows),
            validation=validation,
            input_paths=(metadata_path, coverage_path),
            dpi=dpi,
        )
    finally:
        try:
            import matplotlib.pyplot as plt

            plt.close(figure)
        except ImportError:
            pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--coverage", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="jcvi_bidirectional_coverage")
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        bundle = render_bundle(
            metadata_path=args.metadata,
            coverage_path=args.coverage,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (OSError, JCVIPlotError, ValueError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Published bidirectional JCVI coverage bundle: {bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
