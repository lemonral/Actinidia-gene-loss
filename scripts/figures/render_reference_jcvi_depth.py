#!/usr/bin/env python3
"""Render bidirectional JCVI collinear-block gene-depth coverage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile


SCRIPT_VERSION = "1.1.0"
INPUT_STATUS = "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_GENE_DEPTH_SUMMARY"
OUTPUT_STATUS = "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_GENE_DEPTH_FIGURE"
REQUIRED_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "reference_bed_genes",
    "reference_nonzero_depth_genes",
    "reference_gene_depth_coverage_percent",
    "target_bed_genes",
    "target_nonzero_depth_genes",
    "target_gene_depth_coverage_percent",
    "metric_definition",
)
CHROMOSOME_REQUIRED_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "reference_chromosome",
    "reference_chromosome_genes",
    "reference_nonzero_depth_genes",
    "reference_gene_depth_coverage_percent",
)


class FigureError(RuntimeError):
    """Raised when the JCVI gene-depth figure cannot be validated."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise FigureError(f"{label} is missing or empty: {resolved}")
    return resolved


def read_summary(path: Path) -> list[dict[str, str]]:
    with require_file(path, "gene-depth summary").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise FigureError("gene-depth summary has no header")
        missing = sorted(set(REQUIRED_COLUMNS).difference(reader.fieldnames))
        if missing:
            raise FigureError(f"gene-depth summary is missing columns: {missing}")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    if len(rows) != 23:
        raise FigureError(f"expected 23 gene-depth rows, found {len(rows)}")
    rows.sort(key=lambda row: int(row["plot_order"]))
    if [int(row["plot_order"]) for row in rows] != list(range(1, 24)):
        raise FigureError("plot_order must contain each integer from 1 through 23")
    units = [row["assembly_unit_id"] for row in rows]
    if len(units) != len(set(units)):
        raise FigureError("assembly_unit_id values are not unique")
    return rows


def read_chromosome_summary(path: Path) -> list[dict[str, str]]:
    with require_file(path, "chromosome gene-depth summary").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise FigureError("chromosome gene-depth summary has no header")
        missing = sorted(set(CHROMOSOME_REQUIRED_COLUMNS).difference(reader.fieldnames))
        if missing:
            raise FigureError(
                f"chromosome gene-depth summary is missing columns: {missing}"
            )
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    if len(rows) != 23 * 24:
        raise FigureError(
            f"expected 552 chromosome gene-depth rows, found {len(rows)}"
        )
    return rows


def numeric(row: dict[str, str], column: str) -> float:
    try:
        value = float(row[column])
    except (KeyError, ValueError) as error:
        raise FigureError(f"invalid numeric value in {column}") from error
    if not math.isfinite(value):
        raise FigureError(f"non-finite numeric value in {column}")
    return value


def integer(row: dict[str, str], column: str) -> int:
    value = numeric(row, column)
    if value < 0 or value != int(value):
        raise FigureError(f"{column} must contain nonnegative integers")
    return int(value)


def display_label(species: str, suffix: str) -> str:
    parts = species.split()
    if len(parts) < 2:
        raise FigureError(f"invalid biological species label: {species!r}")
    if len(parts) >= 3 and parts[1] in {"×", "x"}:
        species_text = rf"$\mathit{{{parts[0][0]}.\ {' '.join(parts[2:])}}}$"
    else:
        species_text = rf"$\mathit{{{parts[0][0]}.\ {' '.join(parts[1:])}}}$"
    if suffix.strip().casefold() in {
        "",
        "unphased",
        "actinidiabase v1",
        "unresolved_polyploid_unit",
    }:
        return species_text
    suffix_math = suffix.replace(" ", r"\ ")
    return species_text + rf" $\mathrm{{{suffix_math}}}$"


def validate_rows(rows: list[dict[str, str]]) -> tuple[list[float], list[float]]:
    definitions = {row["metric_definition"] for row in rows}
    if len(definitions) != 1 or not next(iter(definitions)):
        raise FigureError("all rows must share one nonempty metric definition")
    reference_values: list[float] = []
    target_values: list[float] = []
    for row in rows:
        for prefix in ("reference", "target"):
            covered = integer(row, f"{prefix}_nonzero_depth_genes")
            total = integer(row, f"{prefix}_bed_genes")
            percent = numeric(row, f"{prefix}_gene_depth_coverage_percent")
            if total == 0 or covered > total:
                raise FigureError(f"{row['assembly_unit_id']}: invalid {prefix} counts")
            expected = covered * 100.0 / total
            if abs(percent - expected) > 1e-6:
                raise FigureError(
                    f"{row['assembly_unit_id']}: {prefix} percentage does not "
                    "reconcile to covered/total"
                )
            if not 0 <= percent <= 100:
                raise FigureError(f"{prefix} percentage lies outside [0,100]")
            (reference_values if prefix == "reference" else target_values).append(
                percent
            )
    return reference_values, target_values


def render(
    summary_path: Path,
    chromosome_path: Path,
    validation_path: Path,
    output_dir: Path,
    dpi: int,
) -> None:
    summary_path = require_file(summary_path, "gene-depth summary")
    chromosome_path = require_file(
        chromosome_path, "chromosome gene-depth summary"
    )
    validation_path = require_file(validation_path, "summary validation")
    rows = read_summary(summary_path)
    chromosome_rows = read_chromosome_summary(chromosome_path)
    input_validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if (
        input_validation.get("status") not in {INPUT_STATUS, OUTPUT_STATUS}
        or input_validation.get("unit_count") != 23
        or input_validation.get("reference_chromosome_count") != 24
        or input_validation.get("chromosome_grid_rows") != 23 * 24
    ):
        raise FigureError(
            "gene-depth summary validation did not close at 23 units and 24 "
            "reference chromosomes"
        )
    reference_values, target_values = validate_rows(rows)
    units = [row["assembly_unit_id"] for row in rows]
    chromosomes = []
    for chromosome_row in chromosome_rows:
        chromosome_id = chromosome_row["reference_chromosome"]
        if chromosome_id not in chromosomes:
            chromosomes.append(chromosome_id)
        covered = integer(chromosome_row, "reference_nonzero_depth_genes")
        total = integer(chromosome_row, "reference_chromosome_genes")
        percentage = numeric(
            chromosome_row, "reference_gene_depth_coverage_percent"
        )
        if total == 0 or covered > total:
            raise FigureError("invalid chromosome gene-depth counts")
        if abs(percentage - covered * 100.0 / total) > 1e-6:
            raise FigureError("chromosome percentage does not reconcile to counts")
    if len(chromosomes) != 24:
        raise FigureError("chromosome summary must contain 24 reference sequences")
    chromosome_lookup = {
        (row["assembly_unit_id"], row["reference_chromosome"]): numeric(
            row, "reference_gene_depth_coverage_percent"
        )
        for row in chromosome_rows
    }
    if set(chromosome_lookup) != {
        (unit, chromosome) for unit in units for chromosome in chromosomes
    }:
        raise FigureError("chromosome heat-map grid is incomplete")
    heat = [
        [chromosome_lookup[(unit, chromosome)] for chromosome in chromosomes]
        for unit in units
    ]

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError as error:
        raise FigureError("matplotlib is required") from error

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    labels = [
        display_label(row["biological_species"], row["haplotype_or_subgenome"])
        for row in rows
    ]
    y = list(range(len(rows)))
    figure = plt.figure(figsize=(12.2, 8.2))
    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=(0.92, 1.38),
        left=0.205,
        right=0.95,
        bottom=0.105,
        top=0.945,
        wspace=0.12,
    )
    axis = figure.add_subplot(grid[0, 0])
    heat_axis = figure.add_subplot(grid[0, 1], sharey=axis)

    for index, (reference_value, target_value) in enumerate(
        zip(reference_values, target_values, strict=True)
    ):
        axis.plot(
            [target_value, reference_value],
            [index, index],
            color="#B8B8B8",
            linewidth=0.9,
            zorder=1,
        )
    axis.scatter(
        reference_values,
        y,
        s=31,
        marker="o",
        color="#3B6A9A",
        edgecolor="white",
        linewidth=0.5,
        label=r"$\mathit{C.\ scandens}$ genes",
        zorder=3,
    )
    axis.scatter(
        target_values,
        y,
        s=32,
        marker="s",
        color="#D88A3D",
        edgecolor="white",
        linewidth=0.5,
        label="Actinidia-unit genes",
        zorder=3,
    )
    axis.set_xlim(96.5, 100.05)
    axis.set_xticks([97, 98, 99, 100])
    axis.set_ylim(len(rows) - 0.5, -0.5)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Genes covered by JCVI collinear-block depth (%)")
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.005),
        ncol=2,
        borderaxespad=0.0,
        handletextpad=0.45,
        columnspacing=1.2,
    )
    axis.text(
        -0.22,
        1.015,
        "(a)",
        transform=axis.transAxes,
        fontsize=11,
        fontweight="normal",
        va="bottom",
    )

    heat_values = [value for row in heat for value in row]
    heat_minimum = min(heat_values)
    heat_vmin = max(90.0, float(math.floor(heat_minimum)))
    if heat_vmin >= 99.0:
        heat_vmin = 98.0
    colormap = LinearSegmentedColormap.from_list(
        "jcvi_gene_depth", ["#F6E8C3", "#80B8B0", "#2A5F85"]
    )
    image = heat_axis.imshow(
        heat,
        aspect="auto",
        interpolation="nearest",
        cmap=colormap,
        vmin=heat_vmin,
        vmax=100,
    )
    heat_axis.set_xticks(range(len(chromosomes)), chromosomes, rotation=90)
    heat_axis.tick_params(axis="y", left=False, labelleft=False)
    heat_axis.tick_params(axis="x", length=0)
    heat_axis.set_xlabel(r"$\mathit{C.\ scandens}$ reference assembly sequence")
    for spine in heat_axis.spines.values():
        spine.set_linewidth(0.7)
        spine.set_color("#4D4D4D")
    colorbar = figure.colorbar(
        image, ax=heat_axis, fraction=0.035, pad=0.025, aspect=28
    )
    colorbar.set_label("Reference chromosome gene-depth coverage (%)")
    colorbar.set_ticks([heat_vmin, (heat_vmin + 100.0) / 2.0, 100.0])
    heat_axis.text(
        -0.04,
        1.015,
        "(b)",
        transform=heat_axis.transAxes,
        fontsize=11,
        fontweight="normal",
        va="bottom",
    )

    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FigureError(f"refusing existing output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    basename = "clematoclethra_actinidia_jcvi_gene_depth"
    try:
        png = temporary / f"{basename}.png"
        pdf = temporary / f"{basename}.pdf"
        figure.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
        figure.savefig(pdf, bbox_inches="tight", facecolor="white")
        plt.close(figure)
        shutil.copy2(summary_path, temporary / f"{basename}.summary.tsv")
        shutil.copy2(
            chromosome_path,
            temporary / f"{basename}.chromosome_summary.tsv",
        )
        caption = (
            "Bidirectional JCVI collinear-block gene-depth coverage between "
            "Clematoclethra scandens and 23 Actinidia assembly units. (a) Overall "
            "gene-depth coverage calculated independently for the C. scandens and "
            "Actinidia sides. (b) C. scandens gene-depth coverage resolved across "
            "its 24 reference assembly sequences. For each raw JCVI anchor block, "
            "the minimum and maximum ordered BED-gene indices defined a half-open "
            "interval; overlapping intervals were merged, and genes with non-zero "
            "block depth were divided by all BED genes on the same genome side. "
            "Values represent gene-based collinear-block coverage, not nucleotide "
            "sequence identity or base-pair coverage."
        )
        caption_path = temporary / f"{basename}.caption.txt"
        caption_path.write_text(caption + "\n", encoding="utf-8")
        output_records = []
        for path in (
            png,
            pdf,
            temporary / f"{basename}.summary.tsv",
            temporary / f"{basename}.chromosome_summary.tsv",
            caption_path,
        ):
            output_records.append(
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
        figure_validation = {
            "schema_version": 1,
            "status": OUTPUT_STATUS,
            "renderer": "scripts/figures/render_reference_jcvi_depth.py",
            "renderer_version": SCRIPT_VERSION,
            "unit_count": 23,
            "reference_chromosome_count": 24,
            "chromosome_grid_rows": 23 * 24,
            "minimum_reference_gene_depth_coverage_percent": min(reference_values),
            "maximum_reference_gene_depth_coverage_percent": max(reference_values),
            "minimum_target_gene_depth_coverage_percent": min(target_values),
            "maximum_target_gene_depth_coverage_percent": max(target_values),
            "inputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in (summary_path, chromosome_path, validation_path)
            ],
            "outputs": output_records,
            "checks": {
                "twenty_three_units_complete": "pass",
                "complete_23_by_24_chromosome_heatmap": "pass",
                "coverage_reconciles_to_nonzero_depth_over_bed_genes": "pass",
                "gene_depth_labeled_not_sequence_identity": "pass",
                "italic_binomials_upright_suffixes": "pass",
                "subplot_titles_absent": "pass",
            },
        }
        (temporary / f"{basename}.validation.json").write_text(
            json.dumps(figure_validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output_dir)
    except Exception:
        try:
            plt.close(figure)
        except Exception:
            pass
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--chromosome-summary", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args()
    try:
        render(
            args.summary,
            args.chromosome_summary,
            args.validation,
            args.output_dir,
            args.dpi,
        )
        print(f"PASS: {args.output_dir}")
        return 0
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError, FigureError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
