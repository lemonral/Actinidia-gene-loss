#!/usr/bin/env python3
"""Render current OrthoFinder core/species-specific and gene-composition panels.

The four gene-composition classes are mutually exclusive:
  1. genes in orthogroups found only in the focal taxon;
  2. genes present as one focal copy in an orthogroup shared by >=2 taxa;
  3. genes present as >=2 focal copies in an orthogroup shared by >=2 taxa;
  4. genes not assigned to an orthogroup by OrthoFinder.

This avoids the dimensional overlap in the historical plot, where a count of
species-specific families was stacked together with counts of genes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


SCRIPT_VERSION = "1.0.0"

TAXON_ORDER = (
    "act_arguta_C",
    "act_chinensis_A",
    "act_deliciosa_ADM_D",
    "act_eriantha_HAP1",
    "act_hemsleyana",
    "act_latifolia",
    "act_longicarpa",
    "act_macrosperma",
    "act_polygama",
    "act_reticulata",
    "act_rufa_ActinidiaBase_v1",
    "act_zhejiangensis_A",
    "act_zhejiangensis_B",
    "clem_scandens",
    "coffea_arabica_E",
    "rhodo_simsii",
    "vitis_vinifera",
)

DISPLAY = {
    "act_arguta_C": (r"$\it{A.\ arguta}$", "A. arguta"),
    "act_chinensis_A": (r"$\it{A.\ chinensis}$", "A. chinensis"),
    "act_deliciosa_ADM_D": (r"$\it{A.\ deliciosa}$", "A. deliciosa"),
    "act_eriantha_HAP1": (r"$\it{A.\ eriantha}$", "A. eriantha"),
    "act_hemsleyana": (r"$\it{A.\ hemsleyana}$", "A. hemsleyana"),
    "act_latifolia": (r"$\it{A.\ latifolia}$", "A. latifolia"),
    "act_longicarpa": (r"$\it{A.\ longicarpa}$", "A. longicarpa"),
    "act_macrosperma": (r"$\it{A.\ macrosperma}$", "A. macrosperma"),
    "act_polygama": (r"$\it{A.\ polygama}$", "A. polygama"),
    "act_reticulata": (r"$\it{A.\ reticulata}$", "A. reticulata"),
    "act_rufa_ActinidiaBase_v1": (r"$\it{A.\ rufa}$", "A. rufa"),
    "act_zhejiangensis_A": (
        r"$\it{A.\ zhejiangensis}$ A",
        "A. zhejiangensis A",
    ),
    "act_zhejiangensis_B": (
        r"$\it{A.\ zhejiangensis}$ B",
        "A. zhejiangensis B",
    ),
    "clem_scandens": (r"$\it{C.\ scandens}$", "C. scandens"),
    "coffea_arabica_E": (r"$\it{C.\ arabica}$", "C. arabica"),
    "rhodo_simsii": (r"$\it{R.\ simsii}$", "R. simsii"),
    "vitis_vinifera": (r"$\it{V.\ vinifera}$", "V. vinifera"),
}

CATEGORY_COLUMNS = (
    "species_specific_orthogroup_genes",
    "shared_single_copy_orthogroup_genes",
    "shared_multi_copy_orthogroup_genes",
    "unassigned_genes",
)

CATEGORY_LABELS = (
    "Species-specific orthogroup genes",
    "Single-copy genes in shared orthogroups",
    "Multi-copy genes in shared orthogroups",
    "Unassigned genes",
)

CATEGORY_COLORS = ("#4C78A8", "#72B7B2", "#F2CF63", "#E07B62")


class ProfileError(RuntimeError):
    """Raised when the exact OrthoFinder inputs do not reconcile."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_validation(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("status") != "PASS":
        raise ProfileError(f"OrthoFinder validation is not PASS: {data.get('status')!r}")
    return data


def read_gene_counts(path: Path) -> tuple[list[str], list[list[int]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 3 or header[0] != "Orthogroup" or header[-1] != "Total":
            raise ProfileError(f"unexpected GeneCount header: {header!r}")
        taxa = header[1:-1]
        rows: list[list[int]] = []
        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise ProfileError(f"GeneCount row {line_number}: wrong number of fields")
            values = [int(value) for value in row[1:-1]]
            if sum(values) != int(row[-1]):
                raise ProfileError(f"GeneCount row {line_number}: Total does not reconcile")
            rows.append(values)
    if tuple(taxa) != TAXON_ORDER:
        raise ProfileError(
            "current OrthoFinder taxon order/set does not match the frozen 17-taxon design"
        )
    return taxa, rows


def read_statistics(path: Path, taxa: list[str]) -> dict[str, dict[str, int]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != "" or header[1:] != taxa:
            raise ProfileError("Statistics_PerSpecies taxon header does not match GeneCount")
        wanted = {
            "Number of genes": "total_genes",
            "Number of genes in orthogroups": "assigned_genes",
            "Number of unassigned genes": "unassigned_genes",
            "Number of species-specific orthogroups": "reported_species_specific_orthogroups",
            "Number of genes in species-specific orthogroups": (
                "reported_species_specific_orthogroup_genes"
            ),
        }
        result = {taxon: {} for taxon in taxa}
        for row in reader:
            if not row or row[0] not in wanted:
                continue
            if len(row) != len(taxa) + 1:
                raise ProfileError(f"Statistics_PerSpecies row {row[0]!r}: wrong width")
            key = wanted[row[0]]
            for taxon, value in zip(taxa, row[1:]):
                result[taxon][key] = int(value)
    for taxon, values in result.items():
        missing = set(wanted.values()) - set(values)
        if missing:
            raise ProfileError(f"Statistics_PerSpecies {taxon}: missing {sorted(missing)}")
    return result


def summarize(
    taxa: list[str],
    rows: Iterable[list[int]],
    statistics: dict[str, dict[str, int]],
) -> tuple[list[dict[str, int | str]], int]:
    summary = {
        taxon: {
            "species_specific_orthogroups": 0,
            "species_specific_orthogroup_genes": 0,
            "shared_single_copy_orthogroup_genes": 0,
            "shared_multi_copy_orthogroup_genes": 0,
        }
        for taxon in taxa
    }
    all_species_orthogroups = 0
    for values in rows:
        present = sum(value > 0 for value in values)
        if present == len(taxa):
            all_species_orthogroups += 1
        for taxon, value in zip(taxa, values):
            if value == 0:
                continue
            target = summary[taxon]
            if present == 1:
                target["species_specific_orthogroups"] += 1
                target["species_specific_orthogroup_genes"] += value
            elif value == 1:
                target["shared_single_copy_orthogroup_genes"] += 1
            else:
                target["shared_multi_copy_orthogroup_genes"] += value

    result: list[dict[str, int | str]] = []
    for taxon in taxa:
        row = {
            "taxon_id": taxon,
            "display_label": DISPLAY[taxon][1],
            **summary[taxon],
            "unassigned_genes": statistics[taxon]["unassigned_genes"],
            "total_genes": statistics[taxon]["total_genes"],
        }
        if (
            row["species_specific_orthogroups"]
            != statistics[taxon]["reported_species_specific_orthogroups"]
        ):
            raise ProfileError(f"{taxon}: species-specific orthogroup count mismatch")
        if (
            row["species_specific_orthogroup_genes"]
            != statistics[taxon]["reported_species_specific_orthogroup_genes"]
        ):
            raise ProfileError(f"{taxon}: species-specific gene count mismatch")
        assigned = sum(int(row[column]) for column in CATEGORY_COLUMNS[:-1])
        if assigned != statistics[taxon]["assigned_genes"]:
            raise ProfileError(f"{taxon}: assigned-gene categories do not reconcile")
        if sum(int(row[column]) for column in CATEGORY_COLUMNS) != row["total_genes"]:
            raise ProfileError(f"{taxon}: all gene categories do not sum to total genes")
        result.append(row)
    return result, all_species_orthogroups


def configure_matplotlib() -> None:
    import matplotlib as mpl

    mpl.use("Agg")
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )


def _draw_flower(ax, profiles: list[dict[str, int | str]], common: int) -> None:
    from matplotlib.patches import Circle, Ellipse

    n = len(profiles)
    angles = [math.pi / 2 + 2 * math.pi * index / n for index in range(n)]
    petal_colors = ("#DCE6F2", "#DCEDEB", "#F8E8C2", "#F5DCD6")
    for index, (angle, profile) in enumerate(zip(angles, profiles)):
        cx, cy = 1.62 * math.cos(angle), 1.62 * math.sin(angle)
        petal = Ellipse(
            (cx, cy),
            width=1.48,
            height=0.45,
            angle=math.degrees(angle),
            facecolor=petal_colors[index % len(petal_colors)],
            edgecolor="#666666",
            linewidth=0.55,
            zorder=1,
        )
        ax.add_patch(petal)
        count = int(profile["species_specific_orthogroups"])
        ax.text(
            1.88 * math.cos(angle),
            1.88 * math.sin(angle),
            f"{count:,}",
            ha="center",
            va="center",
            fontsize=8.1,
            zorder=3,
        )

        label_radius = 2.55
        x, y = label_radius * math.cos(angle), label_radius * math.sin(angle)
        cosine = math.cos(angle)
        if cosine > 0.22:
            ha = "left"
        elif cosine < -0.22:
            ha = "right"
        else:
            ha = "center"
        sine = math.sin(angle)
        va = "bottom" if sine > 0.35 else "top" if sine < -0.35 else "center"
        ax.text(
            x,
            y,
            DISPLAY[str(profile["taxon_id"])][0],
            ha=ha,
            va=va,
            fontsize=8.2,
        )

    circle = Circle(
        (0, 0),
        radius=0.92,
        facecolor="#F28E2B",
        edgecolor="#333333",
        linewidth=0.85,
        zorder=2,
    )
    ax.add_patch(circle)
    ax.text(
        0,
        0.08,
        f"{common:,}",
        color="white",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        zorder=4,
    )
    ax.text(
        0,
        -0.25,
        "orthogroups present\nin all 17 taxa",
        color="white",
        ha="center",
        va="center",
        fontsize=8.2,
        linespacing=1.2,
        zorder=4,
    )
    ax.set_xlim(-3.25, 3.25)
    ax.set_ylim(-3.25, 3.25)
    ax.set_aspect("equal")
    ax.axis("off")


def _draw_stacked(ax, profiles: list[dict[str, int | str]]) -> None:
    import numpy as np

    ordered = list(reversed(profiles))
    y = np.arange(len(ordered))
    left = np.zeros(len(ordered), dtype=float)
    for column, label, color in zip(
        CATEGORY_COLUMNS, CATEGORY_LABELS, CATEGORY_COLORS
    ):
        values = np.array([int(row[column]) for row in ordered])
        ax.barh(
            y,
            values,
            left=left,
            height=0.72,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            label=label,
        )
        left += values

    ax.set_yticks(y)
    ax.set_yticklabels([DISPLAY[str(row["taxon_id"])][0] for row in ordered])
    ax.set_xlabel("Number of genes")
    ax.set_ylim(-0.7, len(ordered) - 0.3)
    ax.set_xlim(0, max(left) * 1.095)
    ax.xaxis.grid(True, color="#D9D9D9", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for yi, total in zip(y, left):
        ax.text(
            total + max(left) * 0.012,
            yi,
            f"{int(total):,}",
            va="center",
            ha="left",
            fontsize=7.5,
            color="#333333",
        )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.105),
        ncol=2,
        frameon=False,
        columnspacing=1.3,
        handlelength=1.6,
    )


def render(
    profiles: list[dict[str, int | str]],
    common: int,
    outdir: Path,
    prefix: str,
) -> list[Path]:
    configure_matplotlib()
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    figure = plt.figure(figsize=(15.6, 8.2))
    grid = figure.add_gridspec(1, 2, width_ratios=(1.05, 1.0), wspace=0.2)
    ax_a = figure.add_subplot(grid[0, 0])
    ax_b = figure.add_subplot(grid[0, 1])
    _draw_flower(ax_a, profiles, common)
    _draw_stacked(ax_b, profiles)
    ax_a.text(-0.045, 1.025, "(a)", transform=ax_a.transAxes, fontsize=13)
    ax_b.text(-0.11, 1.025, "(b)", transform=ax_b.transAxes, fontsize=13)
    figure.subplots_adjust(left=0.035, right=0.985, top=0.97, bottom=0.15)
    for suffix in ("png", "pdf"):
        path = outdir / f"{prefix}.{suffix}"
        figure.savefig(path, dpi=400 if suffix == "png" else None, bbox_inches="tight")
        outputs.append(path)
    plt.close(figure)

    for panel, drawer, size in (
        ("core_species_specific", lambda ax: _draw_flower(ax, profiles, common), (8.4, 8.0)),
        ("gene_category_composition", lambda ax: _draw_stacked(ax, profiles), (8.3, 8.0)),
    ):
        fig, ax = plt.subplots(figsize=size)
        drawer(ax)
        fig.subplots_adjust(
            left=0.22 if panel == "gene_category_composition" else 0.04,
            right=0.97,
            top=0.98,
            bottom=0.16 if panel == "gene_category_composition" else 0.04,
        )
        for suffix in ("png", "pdf"):
            path = outdir / f"{panel}.{suffix}"
            fig.savefig(path, dpi=400 if suffix == "png" else None, bbox_inches="tight")
            outputs.append(path)
        plt.close(fig)
    return outputs


def write_table(path: Path, profiles: list[dict[str, int | str]]) -> None:
    columns = (
        "taxon_id",
        "display_label",
        "total_genes",
        "species_specific_orthogroups",
        *CATEGORY_COLUMNS,
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(profiles)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-counts", type=Path, required=True)
    parser.add_argument("--statistics", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="orthofinder_species_profiles")
    args = parser.parse_args()

    gene_counts = args.gene_counts.expanduser().resolve()
    statistics_path = args.statistics.expanduser().resolve()
    validation_path = args.validation.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()

    validation = read_validation(validation_path)
    taxa, rows = read_gene_counts(gene_counts)
    statistics = read_statistics(statistics_path, taxa)
    profiles, common = summarize(taxa, rows, statistics)

    if validation.get("species") != len(taxa):
        raise ProfileError("validation species count does not reconcile")
    if validation.get("orthogroups") != len(rows):
        raise ProfileError("validation orthogroup count does not reconcile")
    if validation.get("all_species_orthogroups") != common:
        raise ProfileError("validation all-species orthogroup count does not reconcile")
    if validation.get("genes") != sum(int(row["total_genes"]) for row in profiles):
        raise ProfileError("validation total-gene count does not reconcile")
    if validation.get("unassigned_genes") != sum(
        int(row["unassigned_genes"]) for row in profiles
    ):
        raise ProfileError("validation unassigned-gene count does not reconcile")

    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / "orthofinder_species_profile.tsv"
    write_table(table_path, profiles)
    output_paths = render(profiles, common, outdir, args.prefix)

    caption = (
        "(a) Orthogroups present in all 17 taxa and species-specific orthogroups "
        "from the frozen current OrthoFinder analysis. (b) Mutually exclusive "
        "gene composition per taxon. Species-specific orthogroup genes occur in "
        "orthogroups restricted to one taxon; single-copy and multi-copy genes "
        "occur in orthogroups shared by at least two taxa; unassigned genes were "
        "not placed in an orthogroup by OrthoFinder."
    )
    caption_path = outdir / "caption.txt"
    caption_path.write_text(caption + "\n", encoding="utf-8")

    all_outputs = [table_path, caption_path, *output_paths]
    manifest = {
        "status": "PASS_ORTHOFINDER_SPECIES_PROFILES",
        "schema_version": 1,
        "script_version": SCRIPT_VERSION,
        "taxa": len(taxa),
        "orthogroups": len(rows),
        "all_species_orthogroups": common,
        "total_genes": sum(int(row["total_genes"]) for row in profiles),
        "assigned_genes": sum(
            sum(int(row[column]) for column in CATEGORY_COLUMNS[:-1])
            for row in profiles
        ),
        "unassigned_genes": sum(int(row["unassigned_genes"]) for row in profiles),
        "input_checksums": {
            gene_counts.name: sha256(gene_counts),
            statistics_path.name: sha256(statistics_path),
            validation_path.name: sha256(validation_path),
        },
        "outputs": {path.name: sha256(path) for path in all_outputs},
    }
    (outdir / "validation.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
