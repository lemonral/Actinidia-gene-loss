#!/usr/bin/env python3
"""Render JCVI evidence supporting C. scandens as the loss reference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile


SCRIPT_VERSION = "1.1.0"


class FigureError(RuntimeError):
    """Raised when the publication figure cannot be validated."""


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with require_file(path, "input table").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    if not rows:
        raise FigureError(f"input table has no rows: {path}")
    return rows


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
        "unresolved polyploid unit",
        "unresolved_polyploid_unit",
    }:
        return species_text
    suffix_math = suffix.replace(" ", r"\ ")
    return species_text + rf" $\mathrm{{{suffix_math}}}$"


def number(row: dict[str, str], column: str) -> float:
    try:
        value = float(row[column])
    except (KeyError, ValueError) as error:
        raise FigureError(f"invalid numeric value in {column}") from error
    if not math.isfinite(value):
        raise FigureError(f"non-finite numeric value in {column}")
    return value


def render(
    summary_path: Path,
    chromosome_path: Path,
    validation_path: Path,
    output_dir: Path,
    dpi: int,
) -> None:
    summary_path = require_file(summary_path, "summary table")
    chromosome_path = require_file(chromosome_path, "chromosome table")
    validation_path = require_file(validation_path, "summary validation")
    summary = read_tsv(summary_path)
    chromosome = read_tsv(chromosome_path)
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if (
        validation.get("status")
        != "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_REFERENCE_SUPPORT"
        or validation.get("unit_count") != 23
        or validation.get("reference_chromosome_count") != 24
    ):
        raise FigureError("input validation did not close at 23 units and 24 reference sequences")
    summary.sort(key=lambda row: int(row["plot_order"]))
    units = [row["assembly_unit_id"] for row in summary]
    if len(units) != 23 or len(set(units)) != 23:
        raise FigureError("summary must contain 23 unique units")
    chromosomes = sorted({row["reference_chromosome"] for row in chromosome})
    if len(chromosomes) != 24:
        raise FigureError("chromosome table must contain 24 reference sequences")
    chromosome_lookup = {
        (row["assembly_unit_id"], row["reference_chromosome"]): number(
            row, "anchored_reference_gene_percent"
        )
        for row in chromosome
    }
    if len(chromosome_lookup) != 23 * 24:
        raise FigureError("chromosome heat-map grid is incomplete")

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
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    y = list(range(len(summary)))
    reference_coverage = [
        number(row, "reference_anchored_gene_percent") for row in summary
    ]
    target_coverage = [number(row, "target_anchored_gene_percent") for row in summary]
    labels = [
        display_label(row["biological_species"], row["haplotype_or_subgenome"])
        for row in summary
    ]
    heat = [
        [chromosome_lookup[(unit, chromosome_id)] for chromosome_id in chromosomes]
        for unit in units
    ]

    figure = plt.figure(figsize=(12.2, 8.2))
    grid = figure.add_gridspec(
        1, 2, width_ratios=(0.92, 1.38), left=0.205, right=0.95,
        bottom=0.105, top=0.945, wspace=0.12
    )
    coverage_axis = figure.add_subplot(grid[0, 0])
    heat_axis = figure.add_subplot(grid[0, 1], sharey=coverage_axis)

    for index, (reference_value, target_value) in enumerate(
        zip(reference_coverage, target_coverage, strict=True)
    ):
        coverage_axis.plot(
            [target_value, reference_value],
            [index, index],
            color="#B7B7B7",
            linewidth=0.8,
            zorder=1,
        )
    coverage_axis.scatter(
        reference_coverage, y, s=27, marker="o", color="#3B6A9A",
        edgecolor="white", linewidth=0.45, label=r"$\mathit{C.\ scandens}$ genes",
        zorder=3,
    )
    coverage_axis.scatter(
        target_coverage, y, s=28, marker="s", color="#D88A3D",
        edgecolor="white", linewidth=0.45, label="Actinidia genome genes",
        zorder=3,
    )
    coverage_axis.set_xlim(55, 87)
    coverage_axis.set_xticks([55, 60, 65, 70, 75, 80, 85])
    coverage_axis.set_ylim(len(summary) - 0.5, -0.5)
    coverage_axis.set_yticks(y, labels)
    coverage_axis.set_xlabel("Genes in JCVI collinear anchors (%)")
    coverage_axis.grid(axis="x", color="#D9D9D9", linewidth=0.55)
    coverage_axis.set_axisbelow(True)
    coverage_axis.spines[["top", "right"]].set_visible(False)
    coverage_axis.legend(
        frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.005),
        ncol=1, borderaxespad=0.0, handletextpad=0.5
    )
    coverage_axis.text(
        -0.22, 1.015, "(a)", transform=coverage_axis.transAxes,
        fontsize=11, fontweight="normal", va="bottom"
    )

    colormap = LinearSegmentedColormap.from_list(
        "jcvi_reference_coverage", ["#F6E8C3", "#80B8B0", "#2A5F85"]
    )
    image = heat_axis.imshow(
        heat, aspect="auto", interpolation="nearest", cmap=colormap,
        vmin=55, vmax=90
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
    colorbar.set_label(r"Anchored $\mathit{C.\ scandens}$ genes (%)")
    colorbar.set_ticks([55, 60, 70, 80, 90])
    heat_axis.text(
        -0.04, 1.015, "(b)", transform=heat_axis.transAxes,
        fontsize=11, fontweight="normal", va="bottom"
    )

    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FigureError(f"refusing existing output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    basename = "clematoclethra_actinidia_jcvi_support"
    png = temporary / f"{basename}.png"
    pdf = temporary / f"{basename}.pdf"
    figure.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    shutil.copy2(summary_path, temporary / f"{basename}.summary.tsv")
    shutil.copy2(chromosome_path, temporary / f"{basename}.chromosome_summary.tsv")
    caption = (
        "JCVI collinearity between Clematoclethra scandens and the 23 Actinidia "
        "assembly units. (a) Percentages of C. scandens reference genes and genes "
        "in each Actinidia unit represented in raw JCVI collinear anchors. "
        "(b) Percentages of genes on each of the 24 C. scandens reference assembly "
        "sequences represented in anchors with each Actinidia unit. Coverage is "
        "gene based and does not represent nucleotide sequence identity."
    )
    (temporary / f"{basename}.caption.txt").write_text(caption + "\n", encoding="utf-8")
    output_records = []
    for path in (
        png,
        pdf,
        temporary / f"{basename}.summary.tsv",
        temporary / f"{basename}.chromosome_summary.tsv",
        temporary / f"{basename}.caption.txt",
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
        "status": "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_FIGURE",
        "renderer": "scripts/figures/render_reference_jcvi_support.py",
        "renderer_version": SCRIPT_VERSION,
        "unit_count": 23,
        "reference_sequence_count": 24,
        "minimum_reference_anchored_gene_percent": min(reference_coverage),
        "maximum_reference_anchored_gene_percent": max(reference_coverage),
        "minimum_target_anchored_gene_percent": min(target_coverage),
        "maximum_target_anchored_gene_percent": max(target_coverage),
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
            "complete_23_by_24_heatmap": "pass",
            "gene_based_coverage_labeled_not_identity": "pass",
            "italic_binomials_upright_suffixes": "pass",
            "subplot_titles_absent": "pass",
        },
    }
    (temporary / f"{basename}.validation.json").write_text(
        json.dumps(figure_validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, output_dir)


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
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
