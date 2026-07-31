#!/usr/bin/env python3
"""Render a reconciled assembly-and-annotation QC publication bundle.

The renderer joins four TSV inputs by ``assembly_unit_id``: one metadata
table, basic assembly/annotation statistics, genome BUSCO summaries, and
protein BUSCO summaries.  It fails before plotting unless every identifier
set is identical, every identifier is unique, BUSCO provenance is uniform,
and all numeric values are valid.  Metadata rows are retained in their input
order, including assemblies whose decision status is ``excluded``.

Required metadata columns are ``assembly_unit_id``, ``biological_species``,
``haplotype_or_subgenome``, ``accession``, ``assembly_scope``, and
``decision_status``.  Decision status must be exactly ``current``,
``candidate``, or ``excluded``.  The resulting bundle contains PNG and PDF
figures, the exact plotted TSV, an English caption, JSON validation, and a
checksum manifest.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import TaxonLabelError, format_taxon_label_from_metadata


SCRIPT_VERSION = "1.1.0"
DECISION_STATUSES = ("current", "candidate", "excluded")
METADATA_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "accession",
    "assembly_scope",
    "decision_status",
)
BASIC_COLUMNS = ("assembly_unit_id", "genome_total_bp", "gff_gene_count")
BUSCO_COLUMNS = (
    "assembly_unit_id",
    "busco_version",
    "dataset",
    "dataset_creation_date",
    "mode",
    "C_percent",
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
    "genome_size_bp",
    "genome_size_gb",
    "annotated_gene_count",
    "genome_busco_complete_percent",
    "protein_busco_complete_percent",
    "busco_version",
    "busco_dataset",
    "busco_dataset_creation_date",
    "genome_busco_mode",
    "protein_busco_mode",
)
INTEGER_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)\Z")


class QCPlotError(RuntimeError):
    """Raised when publication inputs cannot be reconciled exactly."""


def read_tsv(
    path: Path, required_columns: Sequence[str], table_name: str
) -> list[dict[str, str]]:
    """Read a nonempty TSV after validating required and unique columns."""

    path = Path(path).expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise QCPlotError(f"cannot open {table_name} {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise QCPlotError(f"{path}: missing TSV header")
        fields = [field.strip() for field in reader.fieldnames]
        if any(not field for field in fields) or len(fields) != len(set(fields)):
            raise QCPlotError(f"{path}: TSV column names must be nonempty and unique")
        missing = sorted(set(required_columns).difference(fields))
        if missing:
            raise QCPlotError(
                f"{path}: missing {table_name} columns: {', '.join(missing)}"
            )
        reader.fieldnames = fields
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, 2):
            if None in raw:
                raise QCPlotError(f"{path}:{line_number}: more values than header columns")
            row = {key: (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            rows.append(row)
    if not rows:
        raise QCPlotError(f"{path}: {table_name} contains no data rows")
    return rows


def index_unique(
    rows: Iterable[dict[str, str]], table_name: str
) -> dict[str, dict[str, str]]:
    """Index rows by a nonempty, unique ``assembly_unit_id``."""

    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        assembly_unit_id = row["assembly_unit_id"]
        if not assembly_unit_id:
            raise QCPlotError(f"{table_name}: empty assembly_unit_id")
        if assembly_unit_id in indexed:
            raise QCPlotError(
                f"{table_name}: duplicate assembly_unit_id {assembly_unit_id!r}"
            )
        indexed[assembly_unit_id] = row
    return indexed


def require_identical_ids(
    metadata_ids: set[str], observed_ids: set[str], table_name: str
) -> None:
    """Require an exact table-to-metadata identifier reconciliation."""

    if observed_ids != metadata_ids:
        missing = sorted(metadata_ids.difference(observed_ids))
        extra = sorted(observed_ids.difference(metadata_ids))
        raise QCPlotError(
            f"{table_name} assembly_unit_id set differs from metadata: "
            f"missing={missing}, extra={extra}"
        )


def parse_integer(value: str, *, field: str, assembly_unit_id: str, positive: bool) -> int:
    """Parse one canonical nonnegative integer."""

    if INTEGER_PATTERN.fullmatch(value) is None:
        raise QCPlotError(
            f"{assembly_unit_id}: {field} must be a nonnegative integer, found {value!r}"
        )
    parsed = int(value)
    if positive and parsed == 0:
        raise QCPlotError(f"{assembly_unit_id}: {field} must be greater than zero")
    return parsed


def parse_percent(value: str, *, field: str, assembly_unit_id: str) -> float:
    """Parse a finite percentage in the closed interval 0--100."""

    try:
        parsed = float(value)
    except ValueError as error:
        raise QCPlotError(
            f"{assembly_unit_id}: {field} must be numeric, found {value!r}"
        ) from error
    if not 0.0 <= parsed <= 100.0:
        raise QCPlotError(
            f"{assembly_unit_id}: {field} must be within [0, 100], found {parsed}"
        )
    return parsed


def prepare_plot_rows(
    metadata_path: Path,
    basic_stats_path: Path,
    genome_busco_path: Path,
    protein_busco_path: Path,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Validate inputs and return exact plotted rows plus validation metadata."""

    metadata_rows = read_tsv(metadata_path, METADATA_COLUMNS, "metadata")
    basic_rows = read_tsv(basic_stats_path, BASIC_COLUMNS, "basic-statistics table")
    genome_rows = read_tsv(genome_busco_path, BUSCO_COLUMNS, "genome BUSCO table")
    protein_rows = read_tsv(protein_busco_path, BUSCO_COLUMNS, "protein BUSCO table")

    metadata = index_unique(metadata_rows, "metadata")
    basic = index_unique(basic_rows, "basic-statistics table")
    genome = index_unique(genome_rows, "genome BUSCO table")
    protein = index_unique(protein_rows, "protein BUSCO table")
    metadata_ids = set(metadata)
    require_identical_ids(metadata_ids, set(basic), "basic-statistics table")
    require_identical_ids(metadata_ids, set(genome), "genome BUSCO table")
    require_identical_ids(metadata_ids, set(protein), "protein BUSCO table")

    signatures = {
        (
            row["busco_version"],
            row["dataset"],
            row["dataset_creation_date"],
        )
        for row in (*genome_rows, *protein_rows)
    }
    if len(signatures) != 1 or any(not value for signature in signatures for value in signature):
        raise QCPlotError(
            "genome and protein BUSCO rows must share one nonempty "
            "version/dataset/creation-date signature"
        )
    busco_version, busco_dataset, busco_date = next(iter(signatures))

    plotted: list[dict[str, object]] = []
    for plot_order, metadata_row in enumerate(metadata_rows, 1):
        assembly_unit_id = metadata_row["assembly_unit_id"]
        status = metadata_row["decision_status"]
        if status not in DECISION_STATUSES:
            raise QCPlotError(
                f"{assembly_unit_id}: decision_status must be one of "
                f"{', '.join(DECISION_STATUSES)}, found {status!r}"
            )
        if not metadata_row["assembly_scope"]:
            raise QCPlotError(f"{assembly_unit_id}: assembly_scope must be nonempty")
        if not metadata_row["accession"]:
            raise QCPlotError(f"{assembly_unit_id}: accession must be nonempty")
        try:
            display_label = format_taxon_label_from_metadata(
                metadata_row,
                suffix_fields=("haplotype_or_subgenome",),
                abbreviate_genus=True,
            )
        except TaxonLabelError as error:
            raise QCPlotError(f"{assembly_unit_id}: invalid taxon metadata: {error}") from error

        genome_bp = parse_integer(
            basic[assembly_unit_id]["genome_total_bp"],
            field="genome_total_bp",
            assembly_unit_id=assembly_unit_id,
            positive=True,
        )
        gene_count = parse_integer(
            basic[assembly_unit_id]["gff_gene_count"],
            field="gff_gene_count",
            assembly_unit_id=assembly_unit_id,
            positive=False,
        )
        genome_complete = parse_percent(
            genome[assembly_unit_id]["C_percent"],
            field="genome BUSCO C_percent",
            assembly_unit_id=assembly_unit_id,
        )
        protein_complete = parse_percent(
            protein[assembly_unit_id]["C_percent"],
            field="protein BUSCO C_percent",
            assembly_unit_id=assembly_unit_id,
        )
        if not genome[assembly_unit_id]["mode"] or not protein[assembly_unit_id]["mode"]:
            raise QCPlotError(f"{assembly_unit_id}: BUSCO mode must be nonempty")

        plotted.append(
            {
                "plot_order": plot_order,
                "assembly_unit_id": assembly_unit_id,
                "biological_species": metadata_row["biological_species"],
                "haplotype_or_subgenome": metadata_row["haplotype_or_subgenome"],
                "accession": metadata_row["accession"],
                "assembly_scope": metadata_row["assembly_scope"],
                "decision_status": status,
                "display_label": display_label,
                "genome_size_bp": genome_bp,
                "genome_size_gb": genome_bp / 1_000_000_000.0,
                "annotated_gene_count": gene_count,
                "genome_busco_complete_percent": genome_complete,
                "protein_busco_complete_percent": protein_complete,
                "busco_version": busco_version,
                "busco_dataset": busco_dataset,
                "busco_dataset_creation_date": busco_date,
                "genome_busco_mode": genome[assembly_unit_id]["mode"],
                "protein_busco_mode": protein[assembly_unit_id]["mode"],
            }
        )

    status_counts = Counter(str(row["decision_status"]) for row in plotted)
    validation: dict[str, object] = {
        "schema_version": "1.0",
        "renderer": "scripts/figures/make_qc_figure.py",
        "renderer_version": SCRIPT_VERSION,
        "status": "pass",
        "assembly_unit_count": len(plotted),
        "excluded_assembly_unit_count": status_counts.get("excluded", 0),
        "decision_status_counts": {
            status: status_counts.get(status, 0) for status in DECISION_STATUSES
        },
        "busco_signature": {
            "version": busco_version,
            "dataset": busco_dataset,
            "dataset_creation_date": busco_date,
        },
        "checks": {
            "metadata_ids_unique": "pass",
            "basic_stats_ids_exactly_match_metadata": "pass",
            "genome_busco_ids_exactly_match_metadata": "pass",
            "protein_busco_ids_exactly_match_metadata": "pass",
            "all_metadata_rows_plotted_including_excluded": "pass",
            "decision_status_domain": "pass",
            "numeric_ranges": "pass",
            "busco_signature_uniform": "pass",
            "taxon_labels_metadata_driven": "pass",
        },
    }
    return plotted, validation


def make_figure(plot_rows: Sequence[Mapping[str, object]]):
    """Create the four-panel QC figure; import Matplotlib only when rendering."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise QCPlotError(
            "Matplotlib is required to render the figure; install the project 'plots' extra"
        ) from error

    status_colors = {
        "current": "#2A9D8F",
        "candidate": "#E9C46A",
        "excluded": "#8D99AE",
    }
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.2,
            "legend.fontsize": 7.5,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    metric_specs = (
        ("genome_size_gb", "Genome size (Gb)", lambda value: f"{value:.2f}"),
        ("annotated_gene_count", "Annotated genes", lambda value: f"{int(value):,}"),
        (
            "genome_busco_complete_percent",
            "Genome BUSCO C (%)",
            lambda value: f"{value:.1f}",
        ),
        (
            "protein_busco_complete_percent",
            "Protein BUSCO C (%)",
            lambda value: f"{value:.1f}",
        ),
    )
    count = len(plot_rows)
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(7.2, 9.2),
        sharey=True,
        constrained_layout=True,
    )
    axes = axes.ravel()
    positions = list(range(count))
    labels = [str(row["display_label"]) for row in plot_rows]
    colors = [status_colors[str(row["decision_status"])] for row in plot_rows]
    for panel, axis, (field, title, formatter) in zip("abcd", axes, metric_specs):
        values = [float(row[field]) for row in plot_rows]
        axis.barh(positions, values, color=colors, edgecolor="none", alpha=0.92)
        axis.set_xlabel(title)
        axis.set_ylim(count - 0.5, -0.5)
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.6, alpha=0.8)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        maximum = max(values) if values else 0.0
        padding = maximum * 0.015 if maximum > 0 else 0.1
        axis.set_xlim(0, maximum * 1.22 if maximum > 0 else 1.0)
        for position, value in zip(positions, values):
            axis.text(value + padding, position, formatter(value), va="center", fontsize=6.8)
        axis.text(
            -0.10,
            1.02,
            f"({panel})",
            transform=axis.transAxes,
            fontsize=10.5,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
    for axis in (axes[0], axes[2]):
        axis.set_yticks(positions, labels=labels)
        axis.tick_params(axis="y", labelleft=True)
    for axis in (axes[1], axes[3]):
        axis.tick_params(axis="y", labelleft=False)
    return figure


def caption_for(plot_rows: Sequence[Mapping[str, object]]) -> str:
    """Return an English, data-derived caption for the QC bundle."""

    first = plot_rows[0]
    return (
        "Assembly and annotation quality control for "
        f"{len(plot_rows)} assembly units, including units later excluded from "
        "gene-loss analysis. Panels report decimal genome size (1 Gb = 10^9 bp), "
        "annotated GFF gene count, genome BUSCO complete percentage, and matched "
        "protein-set BUSCO complete percentage. Latin binomials are italic and "
        "haplotype or subgenome suffixes are upright. "
        "Bar colour denotes current, candidate, or excluded status. Assembly scope "
        "and exact values are retained in the plot-data table. BUSCO used version "
        f"{first['busco_version']} and {first['busco_dataset']} "
        f"(dataset creation date {first['busco_dataset_creation_date']})."
    )


def render_bundle(
    *,
    metadata_path: Path,
    basic_stats_path: Path,
    genome_busco_path: Path,
    protein_busco_path: Path,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> FigureBundle:
    """Prepare, render, and atomically publish one QC bundle."""

    plot_rows, validation = prepare_plot_rows(
        metadata_path,
        basic_stats_path,
        genome_busco_path,
        protein_busco_path,
    )
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
            input_paths=(
                metadata_path,
                basic_stats_path,
                genome_busco_path,
                protein_busco_path,
            ),
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
    parser.add_argument("--basic-stats", required=True, type=Path)
    parser.add_argument("--genome-busco", required=True, type=Path)
    parser.add_argument("--protein-busco", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="assembly_annotation_qc")
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        bundle = render_bundle(
            metadata_path=args.metadata,
            basic_stats_path=args.basic_stats,
            genome_busco_path=args.genome_busco,
            protein_busco_path=args.protein_busco,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (OSError, QCPlotError, ValueError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Published assembly/annotation QC bundle: {bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
