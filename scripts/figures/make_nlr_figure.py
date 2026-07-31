#!/usr/bin/env python3
"""Create one publication bundle for one validated NLR analysis cohort.

Each invocation accepts exactly one unit-summary table and writes exactly one
figure bundle.  Run the script separately for the primary cohort and the
A. rufa sensitivity cohort so their assembly scopes and callable denominators
are never pooled.  Biological-species names are italicized and assembly-unit
suffixes remain upright through the shared metadata-driven label formatter.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label_from_metadata


COHORT_ROLES = frozenset({"primary", "a_rufa_sensitivity"})
PLOT_COLUMNS = (
    "analysis_cohort",
    "cohort_role",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "assembly_scope",
    "taxon_label_mathtext",
    "total_nlr_count",
    "positive_reference_nlr_loss_count",
    "callable_reference_nlr_denominator",
    "positive_reference_nlr_loss_percentage",
    "percentage_status",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--basename", required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--abbreviate-genus", action="store_true")
    return parser.parse_args(argv)


def _parse_integer(value: str, *, context: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{context}: expected an integer, found {value!r}") from exc
    if parsed < 0 or str(parsed) != value:
        raise ValueError(f"{context}: expected a canonical non-negative integer, found {value!r}")
    return parsed


def read_unit_summary(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"unit summary must be a regular non-symlink file: {path}")
    required = {
        "analysis_cohort",
        "cohort_role",
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
        "assembly_scope",
        "total_nlr_count",
        "positive_reference_nlr_loss_count",
        "callable_reference_nlr_denominator",
        "positive_reference_nlr_loss_percentage",
        "percentage_status",
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"unit summary is missing columns: {', '.join(sorted(missing))}")
        raw_rows = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"unit summary:{line_number} has extra tab-delimited fields")
            raw_rows.append((line_number, {field: (raw.get(field) or "").strip() for field in fields}))
    if not raw_rows:
        raise ValueError("unit summary has no data rows")

    cohorts = {row["analysis_cohort"] for _, row in raw_rows}
    roles = {row["cohort_role"] for _, row in raw_rows}
    if "" in cohorts or len(cohorts) != 1:
        raise ValueError(
            "unit summary must contain exactly one non-empty analysis_cohort; "
            "primary and A. rufa sensitivity cohorts require separate bundles"
        )
    if len(roles) != 1 or next(iter(roles)) not in COHORT_ROLES:
        raise ValueError(
            "unit summary must contain exactly one cohort_role: primary or a_rufa_sensitivity"
        )
    cohort = next(iter(cohorts))
    role = next(iter(roles))

    seen_units: set[str] = set()
    rows: list[dict[str, Any]] = []
    undefined_count = 0
    for line_number, row in raw_rows:
        unit = row["assembly_unit_id"]
        if not unit:
            raise ValueError(f"unit summary:{line_number}: empty assembly_unit_id")
        if unit in seen_units:
            raise ValueError(f"unit summary:{line_number}: duplicate assembly_unit_id {unit!r}")
        seen_units.add(unit)
        if not row["biological_species"] or not row["assembly_scope"]:
            raise ValueError(f"unit summary:{line_number}: species and assembly scope must be non-empty")
        total = _parse_integer(row["total_nlr_count"], context=f"unit summary:{line_number}:total_nlr_count")
        positive = _parse_integer(
            row["positive_reference_nlr_loss_count"],
            context=f"unit summary:{line_number}:positive_reference_nlr_loss_count",
        )
        denominator = _parse_integer(
            row["callable_reference_nlr_denominator"],
            context=f"unit summary:{line_number}:callable_reference_nlr_denominator",
        )
        if positive > denominator:
            raise ValueError(
                f"unit summary:{line_number}: positive loss count exceeds callable denominator"
            )
        status = row["percentage_status"]
        if denominator == 0:
            if positive != 0 or row["positive_reference_nlr_loss_percentage"] or status != "undefined_zero_denominator":
                raise ValueError(
                    f"unit summary:{line_number}: zero denominator requires zero calls, an empty percentage, "
                    "and undefined_zero_denominator status"
                )
            percentage: float | None = None
            percentage_output: float | str = ""
            undefined_count += 1
        else:
            if status != "defined":
                raise ValueError(f"unit summary:{line_number}: positive denominator requires defined status")
            try:
                percentage = float(row["positive_reference_nlr_loss_percentage"])
            except ValueError as exc:
                raise ValueError(f"unit summary:{line_number}: invalid loss percentage") from exc
            expected = 100.0 * positive / denominator
            if not math.isfinite(percentage) or not math.isclose(percentage, expected, abs_tol=5e-7):
                raise ValueError(
                    f"unit summary:{line_number}: percentage is inconsistent with positive calls / denominator"
                )
            percentage_output = percentage
        label = format_downstream_taxon_label_from_metadata(
            row,
            suffix_fields=("haplotype_or_subgenome",),
            abbreviate_genus=False,
        )
        rows.append(
            {
                "analysis_cohort": cohort,
                "cohort_role": role,
                "assembly_unit_id": unit,
                "biological_species": row["biological_species"],
                "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                "assembly_scope": row["assembly_scope"],
                "taxon_label_mathtext": label,
                "total_nlr_count": total,
                "positive_reference_nlr_loss_count": positive,
                "callable_reference_nlr_denominator": denominator,
                "positive_reference_nlr_loss_percentage": percentage_output,
                "percentage_status": status,
                "_percentage_numeric": percentage,
            }
        )
    validation = {
        "status": "pass",
        "analysis_cohort": cohort,
        "cohort_role": role,
        "assembly_unit_count": len(rows),
        "undefined_percentage_unit_count": undefined_count,
        "checks": {
            "single_cohort_bundle": "pass",
            "unique_assembly_unit_rows": "pass",
            "counts_and_denominators": "pass",
            "percentage_reconciliation": "pass",
            "italic_species_upright_suffix_labels": "pass",
        },
    }
    return rows, validation


def build_figure(rows: Sequence[Mapping[str, Any]], validation: Mapping[str, Any]):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised when optional dependency is absent
        raise RuntimeError(
            "Matplotlib is required for NLR figures; install the project 'plots' optional dependency"
        ) from exc

    labels = [str(row["taxon_label_mathtext"]) for row in rows]
    total = [int(row["total_nlr_count"]) for row in rows]
    positive = [int(row["positive_reference_nlr_loss_count"]) for row in rows]
    percentages = [row["_percentage_numeric"] for row in rows]
    y = list(range(len(rows)))
    height = max(4.0, 0.48 * len(rows) + 1.8)
    figure, (counts_axis, percentage_axis) = plt.subplots(
        1,
        2,
        figsize=(12.5, height),
        gridspec_kw={"width_ratios": (1.65, 1.0), "wspace": 0.08},
        sharey=True,
    )
    offset = 0.19
    counts_axis.barh(
        [value - offset for value in y],
        total,
        height=0.36,
        color="#4C78A8",
        label="Complete NLR repertoire",
    )
    counts_axis.barh(
        [value + offset for value in y],
        positive,
        height=0.36,
        color="#E45756",
        label="Positive reference-NLR loss calls",
    )
    counts_axis.set_yticks(y, labels=labels)
    counts_axis.invert_yaxis()
    counts_axis.set_xlabel("Number of NLR loci or positive loss calls")
    counts_axis.set_ylabel("Assembly unit")
    counts_axis.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    counts_axis.set_axisbelow(True)
    counts_axis.legend(frameon=False, loc="lower right")
    counts_axis.set_title("(a) Repertoire and positive loss calls", loc="left")

    if len(y) != len(percentages):
        raise ValueError("NLR plot positions and percentages have different lengths")
    defined_y = [position for position, value in zip(y, percentages) if value is not None]
    defined_values = [float(value) for value in percentages if value is not None]
    percentage_axis.barh(defined_y, defined_values, height=0.55, color="#72B7B2")
    for position, value in zip(y, percentages):
        if value is None:
            percentage_axis.text(0, position, "NA", va="center", ha="left", color="#555555")
        else:
            percentage_axis.text(
                float(value),
                position,
                f" {float(value):.1f}%",
                va="center",
                ha="left",
                fontsize=8,
            )
    percentage_axis.tick_params(axis="y", labelleft=False)
    percentage_axis.set_xlabel("Positive loss calls / callable reference NLRs (%)")
    percentage_axis.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    percentage_axis.set_axisbelow(True)
    percentage_axis.set_title("(b) Positive-loss percentage", loc="left")
    role_title = (
        "Primary cohort"
        if validation["cohort_role"] == "primary"
        else r"$\mathit{A.\ rufa}$ sensitivity cohort"
    )
    figure.suptitle(f"NLR repertoire and reference-NLR loss: {role_title}", y=1.005)
    return figure


def make_bundle(args: argparse.Namespace) -> FigureBundle:
    rows, validation = read_unit_summary(args.unit_summary)
    if args.abbreviate_genus:
        for row in rows:
            row["taxon_label_mathtext"] = format_downstream_taxon_label_from_metadata(
                row,
                suffix_fields=("haplotype_or_subgenome",),
                abbreviate_genus=True,
            )
    figure = build_figure(rows, validation)
    try:
        plot_rows = [
            {column: row[column] for column in PLOT_COLUMNS}
            for row in rows
        ]
        role_phrase = (
            "primary"
            if validation["cohort_role"] == "primary"
            else "A. rufa sensitivity"
        )
        caption = (
            f"NLR repertoire and positive reference-NLR loss calls for the {role_phrase} cohort "
            f"({validation['analysis_cohort']}). Panel (a) shows the complete NLR repertoire and "
            "the positive reference-NLR loss-call count for each assembly unit. Panel (b) shows "
            "positive calls as a percentage of that unit's explicitly callable reference-NLR "
            "denominator; NA denotes a zero denominator. Biological-species names are italicized "
            "and haplotype or subgenome suffixes are upright. This bundle contains one cohort only; "
            "primary and A. rufa sensitivity denominators are never pooled."
        )
        return write_figure_bundle(
            figure=figure,
            output_dir=args.output_dir,
            basename=args.basename,
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=[args.unit_summary],
            dpi=args.dpi,
        )
    finally:
        try:
            import matplotlib.pyplot as plt

            plt.close(figure)
        except ImportError:  # pragma: no cover
            pass


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        bundle = make_bundle(args)
    except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Published NLR figure bundle: {bundle.directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
