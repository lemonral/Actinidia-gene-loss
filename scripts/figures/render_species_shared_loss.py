#!/usr/bin/env python3
"""Render biological-species shared/non-shared loss evidence.

The input directory must be one complete ``aggregate_species_loss.py`` output.
Technical haplotypes and subgenomes are used by that upstream aggregation but
are deliberately absent from this species-level figure.  Every bar is split
into mutually exclusive species-gene evidence states, and the denominator is
the complete reference-gene universe.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


CATEGORIES = (
    ("shared_positive_complete", "Shared positive-complete", "#0072B2"),
    ("non_shared_positive_complete", "Non-shared positive-complete", "#D55E00"),
    ("positive_partial", "Partial positive evidence", "#E69F00"),
    ("uncertain", "Uncertain", "#CC79A7"),
    ("confidently_not_positive", "Confidently not positive", "#B8B8B8"),
)

PLOT_COLUMNS = (
    "biological_species",
    "display_label",
    "category",
    "category_label",
    "gene_count",
    "reference_gene_denominator",
    "percentage_of_reference_genes",
)


class SpeciesFigureError(ValueError):
    """Raised when aggregate outputs cannot support an unambiguous figure."""


def _read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        raise SpeciesFigureError(f"required aggregate output is missing: {path.name}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or any(not field for field in fields):
            raise SpeciesFigureError(f"{path.name}: missing or invalid TSV header")
        if len(fields) != len(set(fields)):
            raise SpeciesFigureError(f"{path.name}: duplicate TSV column name")
        rows = list(reader)
    if not rows:
        raise SpeciesFigureError(f"{path.name}: no data rows")
    return rows, fields


def _require_columns(path: Path, fields: Iterable[str], required: Iterable[str]) -> None:
    field_set = set(fields)
    missing = sorted(set(required).difference(field_set))
    if missing:
        raise SpeciesFigureError(
            f"{path.name}: missing required column(s): {', '.join(missing)}"
        )


def _parse_nonnegative_integer(value: str, *, context: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise SpeciesFigureError(f"{context}: expected an integer, found {value!r}") from exc
    if parsed < 0:
        raise SpeciesFigureError(f"{context}: expected a nonnegative integer")
    return parsed


def _parse_boolean(value: str, *, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise SpeciesFigureError(
            f"{context}: expected exactly true or false, found {value!r}"
        )
    return normalized == "true"


def prepare_species_plot(
    aggregate_dir: str | Path,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    """Validate aggregate outputs and return the exact rows to be plotted."""

    source = Path(aggregate_dir)
    matrix_path = source / "species_gene_matrix.tsv"
    prevalence_path = source / "species_prevalence.tsv"
    summary_path = source / "species_loss_summary.json"

    matrix, matrix_fields = _read_tsv(matrix_path)
    prevalence, prevalence_fields = _read_tsv(prevalence_path)
    _require_columns(
        matrix_path,
        matrix_fields,
        (
            "reference_gene_id",
            "biological_species",
            "species_gene_status",
            "species_positive_by_rule",
        ),
    )
    _require_columns(
        prevalence_path,
        prevalence_fields,
        (
            "reference_gene_id",
            "biological_species_count",
            "positive_complete_species_count",
            "positive_partial_species_count",
            "not_positive_species_count",
            "uncertain_species_count",
            "shared_positive_complete",
        ),
    )
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise SpeciesFigureError(
            f"required aggregate output is missing: {summary_path.name}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise SpeciesFigureError(f"{summary_path.name}: invalid JSON: {exc}") from exc
    legacy_summary = (
        isinstance(summary, Mapping)
        and "schema_version" not in summary
        and summary.get("status") == "complete"
    )
    current_checks = summary.get("checks") if isinstance(summary, Mapping) else None
    current_summary = (
        isinstance(summary, Mapping)
        and summary.get("schema_version") == "2.0"
        and summary.get("status") == "PASS"
        and isinstance(current_checks, Mapping)
        and bool(current_checks)
        and all(value is True for value in current_checks.values())
    )
    if not (legacy_summary or current_summary):
        raise SpeciesFigureError(
            f"{summary_path.name}: require legacy complete or schema-2.0 PASS with all checks true"
        )

    allowed_statuses = {
        "positive_complete",
        "positive_partial",
        "not_positive",
        "uncertain",
    }
    by_gene: dict[str, dict[str, str]] = defaultdict(dict)
    species: set[str] = set()
    for line_number, row in enumerate(matrix, start=2):
        gene = row["reference_gene_id"].strip()
        taxon = row["biological_species"].strip()
        status = row["species_gene_status"].strip()
        if not gene or not taxon:
            raise SpeciesFigureError(
                f"{matrix_path.name}:{line_number}: gene and biological species are required"
            )
        if status not in allowed_statuses:
            raise SpeciesFigureError(
                f"{matrix_path.name}:{line_number}: unsupported species_gene_status {status!r}"
            )
        _parse_boolean(
            row["species_positive_by_rule"],
            context=f"{matrix_path.name}:{line_number}:species_positive_by_rule",
        )
        if taxon in by_gene[gene]:
            raise SpeciesFigureError(
                f"{matrix_path.name}:{line_number}: duplicate species-gene row for {taxon!r}/{gene!r}"
            )
        by_gene[gene][taxon] = status
        species.add(taxon)

    ordered_species = sorted(species)
    ordered_genes = sorted(by_gene)
    expected_species = set(ordered_species)
    for gene in ordered_genes:
        observed = set(by_gene[gene])
        if observed != expected_species:
            missing = ", ".join(sorted(expected_species.difference(observed)))
            extra = ", ".join(sorted(observed.difference(expected_species)))
            raise SpeciesFigureError(
                f"{matrix_path.name}: incomplete biological-species grid for {gene!r}; "
                f"missing={missing or 'none'}, extra={extra or 'none'}"
            )

    prevalence_by_gene: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(prevalence, start=2):
        gene = row["reference_gene_id"].strip()
        if not gene:
            raise SpeciesFigureError(
                f"{prevalence_path.name}:{line_number}: empty reference_gene_id"
            )
        if gene in prevalence_by_gene:
            raise SpeciesFigureError(
                f"{prevalence_path.name}:{line_number}: duplicate reference gene {gene!r}"
            )
        prevalence_by_gene[gene] = row
    if set(prevalence_by_gene) != set(ordered_genes):
        raise SpeciesFigureError(
            f"{prevalence_path.name}: reference-gene universe does not match species_gene_matrix.tsv"
        )

    shared_by_gene: dict[str, bool] = {}
    for gene in ordered_genes:
        status_counts = Counter(by_gene[gene].values())
        row = prevalence_by_gene[gene]
        reported_species_count = _parse_nonnegative_integer(
            row["biological_species_count"],
            context=f"{prevalence_path.name}:{gene}:biological_species_count",
        )
        if reported_species_count != len(ordered_species):
            raise SpeciesFigureError(
                f"{prevalence_path.name}:{gene}: biological_species_count is "
                f"{reported_species_count}; expected {len(ordered_species)}"
            )
        count_columns = {
            "positive_complete_species_count": "positive_complete",
            "positive_partial_species_count": "positive_partial",
            "not_positive_species_count": "not_positive",
            "uncertain_species_count": "uncertain",
        }
        for column, status in count_columns.items():
            reported = _parse_nonnegative_integer(
                row[column], context=f"{prevalence_path.name}:{gene}:{column}"
            )
            if reported != status_counts[status]:
                raise SpeciesFigureError(
                    f"{prevalence_path.name}:{gene}: {column}={reported} disagrees with "
                    f"species_gene_matrix.tsv ({status_counts[status]})"
                )
        expected_shared = status_counts["positive_complete"] == len(ordered_species)
        reported_shared = _parse_boolean(
            row["shared_positive_complete"],
            context=f"{prevalence_path.name}:{gene}:shared_positive_complete",
        )
        if reported_shared != expected_shared:
            raise SpeciesFigureError(
                f"{prevalence_path.name}:{gene}: shared_positive_complete disagrees with "
                "the complete species matrix"
            )
        shared_by_gene[gene] = reported_shared

    reference_gene_count = len(ordered_genes)
    shared_gene_count = sum(shared_by_gene.values())
    summary_checks = {
        "biological_species_count": len(ordered_species),
        "reference_gene_count": reference_gene_count,
        "shared_positive_complete_gene_count": shared_gene_count,
    }
    for key, expected in summary_checks.items():
        if key in summary and summary[key] != expected:
            raise SpeciesFigureError(
                f"{summary_path.name}: {key}={summary[key]!r}; expected {expected}"
            )

    category_labels = {key: label for key, label, _ in CATEGORIES}
    counts_by_species: dict[str, Counter[str]] = {
        taxon: Counter() for taxon in ordered_species
    }
    for gene in ordered_genes:
        for taxon, status in by_gene[gene].items():
            if status == "positive_complete":
                category = (
                    "shared_positive_complete"
                    if shared_by_gene[gene]
                    else "non_shared_positive_complete"
                )
            elif status == "positive_partial":
                category = "positive_partial"
            elif status == "uncertain":
                category = "uncertain"
            else:
                category = "confidently_not_positive"
            counts_by_species[taxon][category] += 1

    plot_rows: list[dict[str, object]] = []
    for taxon in ordered_species:
        display_label = format_downstream_taxon_label(taxon, abbreviate_genus=True)
        observed_total = sum(counts_by_species[taxon].values())
        if observed_total != reference_gene_count:
            raise SpeciesFigureError(
                f"internal reconciliation failure for {taxon}: {observed_total} != "
                f"{reference_gene_count}"
            )
        for category, _, _ in CATEGORIES:
            count = counts_by_species[taxon][category]
            plot_rows.append(
                {
                    "biological_species": taxon,
                    "display_label": display_label,
                    "category": category,
                    "category_label": category_labels[category],
                    "gene_count": count,
                    "reference_gene_denominator": reference_gene_count,
                    "percentage_of_reference_genes": 100 * count / reference_gene_count,
                }
            )

    validation: dict[str, object] = {
        "schema_version": "1.0",
        "status": "pass",
        "figure_scope": "biological_species_not_technical_assembly_units",
        "biological_species_count": len(ordered_species),
        "reference_gene_count": reference_gene_count,
        "species_gene_matrix_rows": len(matrix),
        "shared_positive_complete_gene_count": shared_gene_count,
        "category_order": [key for key, _, _ in CATEGORIES],
        "checks": {
            "complete_species_gene_grid": True,
            "prevalence_counts_reconciled": True,
            "shared_flag_recomputed": True,
            "summary_counts_reconciled": True,
            "technical_assembly_ids_used_as_axis_labels": False,
        },
    }
    return plot_rows, validation, [matrix_path, prevalence_path, summary_path]


def build_species_figure(plot_rows: list[dict[str, object]]):
    """Build a horizontal percentage figure from validated plotting rows."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SpeciesFigureError(
            "matplotlib is required to render publication figures; install the plots extra"
        ) from exc

    species = []
    labels: dict[str, str] = {}
    by_species_category: dict[tuple[str, str], dict[str, object]] = {}
    for row in plot_rows:
        taxon = str(row["biological_species"])
        if taxon not in labels:
            species.append(taxon)
            labels[taxon] = str(row["display_label"])
        by_species_category[(taxon, str(row["category"]))] = row

    figure_height = max(3.4, 0.44 * len(species) + 1.8)
    figure, axis = plt.subplots(figsize=(10.8, figure_height), constrained_layout=True)
    positions = list(range(len(species)))
    left = [0.0] * len(species)
    for category, label, color in CATEGORIES:
        widths = [
            float(by_species_category[(taxon, category)]["percentage_of_reference_genes"])
            for taxon in species
        ]
        bars = axis.barh(
            positions,
            widths,
            left=left,
            height=0.72,
            label=label,
            color=color,
            edgecolor="white",
            linewidth=0.35,
        )
        for index, (bar, width) in enumerate(zip(bars, widths)):
            count = int(by_species_category[(species[index], category)]["gene_count"])
            if width >= 5 and count:
                text_color = "white" if category in {
                    "shared_positive_complete",
                    "non_shared_positive_complete",
                } else "black"
                axis.text(
                    left[index] + width / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{count}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=text_color,
                )
        left = [start + width for start, width in zip(left, widths)]

    denominator = int(plot_rows[0]["reference_gene_denominator"])
    axis.set_yticks(positions, [labels[taxon] for taxon in species])
    axis.invert_yaxis()
    axis.set_xlim(0, 100)
    axis.set_xlabel("Percentage of the complete reference-gene universe")
    axis.set_ylabel("Biological species")
    axis.set_title(
        f"Species-level shared and non-shared loss evidence (N = {denominator:,} reference genes)",
        pad=9,
    )
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.55)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.12),
        ncol=len(CATEGORIES),
        frameon=False,
        fontsize=7.5,
    )
    return figure


def publish_species_figure(
    *,
    aggregate_dir: str | Path,
    output_dir: str | Path,
    basename: str = "species_shared_nonshared_loss",
    dpi: int = 300,
) -> FigureBundle:
    plot_rows, validation, input_paths = prepare_species_plot(aggregate_dir)
    figure = build_species_figure(plot_rows)
    caption = (
        "Species-level shared and non-shared gene-loss evidence. Each horizontal bar uses "
        "biological species as the unit and partitions the complete reference-gene universe "
        "into mutually exclusive states. Shared positive-complete means that every included "
        "assembly unit of every included biological species has a positive call. Non-shared "
        "positive-complete is complete within the plotted species but not shared by all species. "
        "Partial positive evidence contains a mixture of positive and non-positive technical "
        "assembly units within a species and is not labelled a confident lineage-specific loss. "
        "Uncertain means that no unit is positive and at least one unit is uncertain or not "
        "callable. Numbers inside sufficiently wide segments are gene counts; segment width is "
        "the percentage of the complete reference-gene denominator. Latin binomials are italic."
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
        except ImportError:  # pragma: no cover - handled during figure construction
            pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggregate-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="species_shared_nonshared_loss")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        bundle = publish_species_figure(
            aggregate_dir=args.aggregate_dir,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (OSError, SpeciesFigureError, ValueError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"figure_bundle\t{bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
