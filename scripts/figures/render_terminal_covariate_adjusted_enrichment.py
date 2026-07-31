#!/usr/bin/env python3
"""Render review figures for covariate-adjusted terminal GO/KEGG enrichment."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


CATEGORIES = (
    ("GO biological process", "GO", "biological_process", 12),
    ("GO molecular function", "GO", "molecular_function", 12),
    ("GO cellular component", "GO", "cellular_component", 10),
    ("KEGG orthology", "KEGG_KO", "", 12),
    ("KEGG pathway", "KEGG_PATHWAY", "", 10),
)


class FigureError(ValueError):
    """Raised when adjusted-enrichment output cannot support the figures."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--kegg-pathway-names", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-terminals", type=int, default=23)
    parser.add_argument(
        "--analysis-level",
        choices=("assembly_unit", "biological_species"),
        default="assembly_unit",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise FigureError(f"{path.name}: invalid or empty table")
    return rows


def read_pathway_names(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    expected = {"term_id", "term_name"}
    if not expected.issubset(rows[0]):
        raise FigureError(f"{path.name}: expected columns {sorted(expected)}")
    names = {
        row["term_id"].strip(): row["term_name"].strip()
        for row in rows
        if row["term_id"].strip() and row["term_name"].strip()
    }
    if len(names) != len(rows):
        raise FigureError(f"{path.name}: duplicate or empty pathway mapping")
    return names


def write_tsv(
    path: Path, fields: list[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def display_label(species: str, suffix: str) -> str:
    if "zhejiangensis parental lineage A" in species:
        return r"$\it{A.\ zhejiangensis}$ A"
    if "zhejiangensis parental lineage B" in species:
        return r"$\it{A.\ zhejiangensis}$ B"
    epithet = species.split()[-1]
    base = rf"$\it{{A.\ {epithet}}}$"
    if suffix and suffix != "ActinidiaBase v1":
        return f"{base} {suffix}"
    return base


def category_of(row: Mapping[str, str]) -> str:
    if row["ontology"] == "GO":
        return {
            "biological_process": "GO biological process",
            "molecular_function": "GO molecular function",
            "cellular_component": "GO cellular component",
        }[row["go_namespace"]]
    if row["ontology"] == "KEGG_KO":
        return "KEGG orthology"
    return "KEGG pathway"


def read_go_graph(
    path: Path,
) -> tuple[dict[str, set[str]], dict[str, str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    aliases: dict[str, str] = {}
    current_id = ""
    current_alt: list[str] = []
    obsolete = False

    def store() -> None:
        if current_id and not obsolete:
            aliases[current_id] = current_id
            for alternate in current_alt:
                aliases[alternate] = current_id

    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        if line == "[Term]":
            store()
            current_id = ""
            current_alt = []
            obsolete = False
        elif line.startswith("id: GO:"):
            current_id = line.split(":", 1)[1].strip()
        elif line.startswith("alt_id: GO:"):
            current_alt.append(line.split(":", 1)[1].strip())
        elif line.startswith("is_a: GO:") and current_id:
            parents[current_id].add(line.split()[1])
        elif line == "is_obsolete: true":
            obsolete = True
    store()
    return parents, aliases


def related_checker(
    parents: Mapping[str, set[str]], aliases: Mapping[str, str]
):
    cache: dict[str, set[str]] = {}

    def ancestors(term: str) -> set[str]:
        term = aliases.get(term, term)
        if term in cache:
            return cache[term]
        result: set[str] = set()
        stack = list(parents.get(term, set()))
        while stack:
            parent = stack.pop()
            parent = aliases.get(parent, parent)
            if parent in result:
                continue
            result.add(parent)
            stack.extend(parents.get(parent, set()))
        cache[term] = result
        return result

    def related(left: str, right: str) -> bool:
        left = aliases.get(left, left)
        right = aliases.get(right, right)
        return left == right or left in ancestors(right) or right in ancestors(left)

    return related


def choose_terms(
    rows: list[dict[str, str]],
    category: str,
    limit: int,
    related,
) -> list[tuple[str, str]]:
    selected_rows = [
        row
        for row in rows
        if category_of(row) == category
        and row["significant_adjusted"] == "true"
        and row["full_fit_converged"] == "true"
    ]
    by_term: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in selected_rows:
        by_term[(row["term_id"], row["term_name"])].append(row)
    ranked = sorted(
        by_term,
        key=lambda key: (
            -len({row["assembly_unit_id"] for row in by_term[key]}),
            min(float(row["q_score_bh"]) for row in by_term[key]),
            -sum(
                math.log2(float(row["adjusted_odds_ratio"]))
                for row in by_term[key]
            )
            / len(by_term[key]),
            key[0],
        ),
    )
    chosen: list[tuple[str, str]] = []
    chosen_units: dict[str, set[str]] = {}
    for term_id, term_name in ranked:
        units = {row["assembly_unit_id"] for row in by_term[(term_id, term_name)]}
        redundant = False
        if category.startswith("GO "):
            for old_id, _ in chosen:
                old_units = chosen_units[old_id]
                union = units | old_units
                jaccard = len(units & old_units) / len(union) if union else 0.0
                if related(term_id, old_id) and jaccard >= 0.65:
                    redundant = True
                    break
        if redundant:
            continue
        chosen.append((term_id, term_name))
        chosen_units[term_id] = units
        if len(chosen) == limit:
            break
    return chosen


def configure_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "font.size": 9,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(figure, directory: Path, basename: str) -> list[Path]:
    png = directory / f"{basename}.png"
    pdf = directory / f"{basename}.pdf"
    figure.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    figure.savefig(pdf, bbox_inches="tight", facecolor="white")
    return [png, pdf]


def render_overview(
    model_rows: list[dict[str, str]],
    all_rows: list[dict[str, str]],
    output_dir: Path,
) -> list[Path]:
    import matplotlib.pyplot as plt
    import numpy as np

    unit_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in model_rows:
        unit = row["assembly_unit_id"]
        if unit not in seen:
            seen.add(unit)
            unit_rows.append(row)
    units = [row["assembly_unit_id"] for row in unit_rows]
    labels = [
        display_label(row["biological_species"], row["haplotype_or_subgenome"])
        for row in unit_rows
    ]
    category_names = [item[0] for item in CATEGORIES]
    adjusted_counts = np.zeros((len(units), len(category_names)), dtype=int)
    unadjusted_counts = np.zeros_like(adjusted_counts)
    unit_index = {unit: index for index, unit in enumerate(units)}
    category_index = {category: index for index, category in enumerate(category_names)}
    for row in all_rows:
        i = unit_index[row["assembly_unit_id"]]
        j = category_index[category_of(row)]
        if (
            row["significant_adjusted"] == "true"
            and row["full_fit_converged"] == "true"
        ):
            adjusted_counts[i, j] += 1
        if (
            float(row["q_hypergeometric_bh"]) <= 0.05
            and float(row["fold_enrichment_unadjusted"]) > 1
        ):
            unadjusted_counts[i, j] += 1
    foreground = np.asarray(
        [int(row["terminal_event_genes"]) for row in unit_rows], dtype=int
    )
    figure_height = max(5.2, 0.30 * len(units) + 1.3)
    figure = plt.figure(figsize=(7.2, figure_height))
    grid = figure.add_gridspec(
        1, 2, width_ratios=(0.94, 1.18), wspace=0.36
    )
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    y = np.arange(len(units))
    axis_a.barh(y, foreground, color="#547AA5", height=0.72)
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Terminal loss-event genes")
    axis_a.spines[["top", "right"]].set_visible(False)
    axis_a.grid(axis="x", color="#D9D9D9", linewidth=0.5, zorder=0)
    axis_a.set_axisbelow(True)
    axis_a.text(-0.20, 1.012, "(a)", transform=axis_a.transAxes, fontsize=10.5)

    image = axis_b.imshow(
        adjusted_counts,
        cmap="YlOrRd",
        aspect="auto",
        interpolation="nearest",
        vmin=0,
    )
    axis_b.set_yticks(y, labels)
    short_categories = ("GO BP", "GO MF", "GO CC", "KEGG KO", "KEGG pathway")
    axis_b.set_xticks(
        np.arange(len(short_categories)),
        short_categories,
        rotation=42,
        ha="right",
        rotation_mode="anchor",
    )
    for i in range(adjusted_counts.shape[0]):
        for j in range(adjusted_counts.shape[1]):
            value = adjusted_counts[i, j]
            color = "white" if value > adjusted_counts.max() * 0.52 else "#202020"
            axis_b.text(
                j,
                i,
                str(value),
                ha="center",
                va="center",
                fontsize=6.7,
                color=color,
            )
    colorbar = figure.colorbar(image, ax=axis_b, fraction=0.045, pad=0.03)
    colorbar.set_label("Adjusted significant terms", fontsize=8)
    axis_b.text(-0.19, 1.012, "(b)", transform=axis_b.transAxes, fontsize=10.5)
    outputs = save_figure(
        figure, output_dir, "terminal_covariate_adjusted_enrichment_overview"
    )
    plt.close(figure)

    comparison_rows = []
    for i, unit in enumerate(units):
        for j, category in enumerate(category_names):
            comparison_rows.append(
                {
                    "assembly_unit_id": unit,
                    "display_label": labels[i].replace("$", "").replace("\\it", ""),
                    "functional_category": category,
                    "unadjusted_significant_terms": int(unadjusted_counts[i, j]),
                    "adjusted_significant_terms": int(adjusted_counts[i, j]),
                }
            )
    write_tsv(
        output_dir / "significant_term_count_comparison.tsv",
        list(comparison_rows[0]),
        comparison_rows,
    )
    return outputs + [output_dir / "significant_term_count_comparison.tsv"]


def render_detail(
    all_rows: list[dict[str, str]],
    unit_order: list[str],
    unit_labels: list[str],
    category: str,
    selected: list[tuple[str, str]],
    output_dir: Path,
) -> tuple[list[Path], list[dict[str, object]]]:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize

    lookup = {
        (row["assembly_unit_id"], row["term_id"]): row
        for row in all_rows
        if category_of(row) == category
        and row["significant_adjusted"] == "true"
        and row["full_fit_converged"] == "true"
    }
    term_labels = [
        "\n".join(textwrap.wrap(name, width=46)) + f"\n[{term_id}]"
        for term_id, name in selected
    ]
    xs: list[int] = []
    ys: list[int] = []
    colors: list[float] = []
    sizes: list[float] = []
    selected_rows: list[dict[str, object]] = []
    for y, (term_id, term_name) in enumerate(selected):
        for x, unit in enumerate(unit_order):
            row = lookup.get((unit, term_id))
            if row is None:
                continue
            q_value = max(float(row["q_score_bh"]), 1e-300)
            odds_ratio = float(row["adjusted_odds_ratio"])
            xs.append(x)
            ys.append(y)
            colors.append(math.log2(odds_ratio))
            sizes.append(18 + 15 * min(6.0, -math.log10(q_value)))
            selected_rows.append(
                {
                    "functional_category": category,
                    "term_id": term_id,
                    "term_name": term_name,
                    "assembly_unit_id": unit,
                    "terminal_loss_count": row["terminal_loss_count"],
                    "background_count": row["background_count"],
                    "adjusted_odds_ratio": row["adjusted_odds_ratio"],
                    "adjusted_ci95_low": row["adjusted_ci95_low"],
                    "adjusted_ci95_high": row["adjusted_ci95_high"],
                    "q_score_bh": row["q_score_bh"],
                }
            )
    height = max(4.4, 0.48 * len(selected) + 1.8)
    figure, axis = plt.subplots(figsize=(7.2, height))
    norm = Normalize(vmin=0, vmax=max(3.0, min(6.0, max(colors, default=3.0))))
    points = axis.scatter(
        xs,
        ys,
        c=colors,
        s=sizes,
        cmap="viridis",
        norm=norm,
        edgecolors="#222222",
        linewidths=0.25,
    )
    axis.set_xticks(range(len(unit_order)), unit_labels, rotation=62, ha="right")
    axis.set_yticks(range(len(selected)), term_labels)
    axis.invert_yaxis()
    axis.set_xlim(-0.7, len(unit_order) - 0.3)
    axis.grid(color="#E2E2E2", linewidth=0.45)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    colorbar = figure.colorbar(points, ax=axis, fraction=0.03, pad=0.02)
    colorbar.set_label(r"Adjusted $\log_2$(odds ratio)", fontsize=8)
    size_handles = []
    size_labels = []
    for q_value in (0.05, 0.001, 1e-6):
        size = 18 + 15 * min(6.0, -math.log10(q_value))
        size_handles.append(
            axis.scatter(
                [],
                [],
                s=size,
                facecolor="#BDBDBD",
                edgecolor="#222222",
                linewidth=0.25,
            )
        )
        size_labels.append(f"q={q_value:g}")
    axis.legend(
        size_handles,
        size_labels,
        title="BH q",
        loc="lower right",
        bbox_to_anchor=(1.0, 1.01),
        ncol=3,
        frameon=False,
        borderaxespad=0,
        handletextpad=0.25,
        columnspacing=0.75,
        scatterpoints=1,
        fontsize=6.8,
        title_fontsize=6.8,
    )
    basename = {
        "GO biological process": "terminal_adjusted_go_bp",
        "GO molecular function": "terminal_adjusted_go_mf",
        "GO cellular component": "terminal_adjusted_go_cc",
        "KEGG orthology": "terminal_adjusted_kegg_ko",
        "KEGG pathway": "terminal_adjusted_kegg_pathway",
    }[category]
    outputs = save_figure(figure, output_dir, basename)
    plt.close(figure)
    return outputs, selected_rows


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    configure_style()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: output directory exists: {args.output_dir}")
    analysis_manifest = args.analysis_dir / "run_manifest.json"
    model_path = args.analysis_dir / "terminal_model_summary.tsv"
    all_path = args.analysis_dir / "terminal_adjusted_enrichment_all_terms.tsv.gz"
    significant_path = (
        args.analysis_dir / "terminal_adjusted_enrichment_significant.tsv"
    )
    for path in (
        analysis_manifest,
        model_path,
        all_path,
        significant_path,
        args.go_obo,
        args.kegg_pathway_names,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing or empty input: {path}")
    manifest = json.loads(analysis_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_TERMINAL_COVARIATE_ADJUSTED_GO_KEGG":
        raise SystemExit("ERROR: adjusted analysis is not PASS")
    model_rows = read_tsv(model_path)
    all_rows = read_tsv(all_path)
    pathway_names = read_pathway_names(args.kegg_pathway_names)
    for row in all_rows:
        if row["ontology"] == "KEGG_PATHWAY":
            row["term_name"] = pathway_names.get(row["term_id"], row["term_name"])
    if (
        len(model_rows) != 3 * args.expected_terminals
        or len(all_rows) != manifest["counts"]["tested_term_rows"]
    ):
        raise SystemExit("ERROR: adjusted tables do not close against manifest")
    if manifest.get("analysis_level", "assembly_unit") != args.analysis_level:
        raise SystemExit("ERROR: requested analysis level differs from manifest")
    parents, aliases = read_go_graph(args.go_obo)
    related = related_checker(parents, aliases)
    selected_by_category: dict[str, list[tuple[str, str]]] = {}
    for category, _, _, limit in CATEGORIES:
        selected = choose_terms(all_rows, category, limit, related)
        if len(selected) < max(5, limit // 2):
            raise SystemExit(f"ERROR: too few representative terms for {category}")
        selected_by_category[category] = selected

    unit_order: list[str] = []
    unit_labels: list[str] = []
    for row in model_rows:
        unit = row["assembly_unit_id"]
        if unit not in unit_order:
            unit_order.append(unit)
            unit_labels.append(
                display_label(
                    row["biological_species"], row["haplotype_or_subgenome"]
                )
            )
    if len(unit_order) != args.expected_terminals:
        raise SystemExit(
            f"ERROR: expected {args.expected_terminals} terminals in model summary"
        )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.", dir=args.output_dir.parent
        )
    )
    try:
        outputs = render_overview(model_rows, all_rows, temporary)
        selected_rows: list[dict[str, object]] = []
        for category, _, _, _ in CATEGORIES:
            category_outputs, category_rows = render_detail(
                all_rows,
                unit_order,
                unit_labels,
                category,
                selected_by_category[category],
                temporary,
            )
            outputs.extend(category_outputs)
            selected_rows.extend(category_rows)
        selected_path = temporary / "representative_adjusted_terms.tsv"
        write_tsv(
            selected_path,
            [
                "functional_category",
                "term_id",
                "term_name",
                "assembly_unit_id",
                "terminal_loss_count",
                "background_count",
                "adjusted_odds_ratio",
                "adjusted_ci95_low",
                "adjusted_ci95_high",
                "q_score_bh",
            ],
            selected_rows,
        )
        outputs.append(selected_path)
        if args.analysis_level == "biological_species":
            foreground_definition = (
                "Biological-species foregrounds contain article-method complete "
                "losses for which all constituent assembly units are decayed/deleted, "
                "assigned to the focal species terminal after excluding ancestral "
                "loss events."
            )
        else:
            foreground_definition = (
                "Assembly-unit foregrounds contain article-method decayed plus "
                "deleted genes assigned to the focal terminal after excluding "
                "ancestral loss events."
            )
        caption = (
            "Covariate-adjusted functional enrichment among terminal loss events. "
            f"{foreground_definition} Logistic score tests adjust for four-tissue "
            "mean TPM and reference CD-HIT family size; q values are BH-adjusted "
            "within terminal and ontology. The overview reports terminal foreground "
            "sizes and stable adjusted significant-term counts. Detail plots show "
            "nonredundant representative terms; point size encodes significance and "
            "colour encodes the complete-model adjusted odds ratio."
        )
        caption_path = temporary / "caption.txt"
        caption_path.write_text(caption + "\n", encoding="utf-8")
        outputs.append(caption_path)
        figure_manifest = {
            "schema_version": "1.0",
            "status": "PASS_TERMINAL_COVARIATE_ADJUSTED_ENRICHMENT_FIGURES",
            "analysis_status": manifest["status"],
            "analysis_level": args.analysis_level,
            "analysis_counts": manifest["counts"],
            "representative_term_counts": {
                category: len(values)
                for category, values in selected_by_category.items()
            },
            "definitions": {
                "adjusted_significance": (
                    "score-test BH q <= 0.05, positive association, and converged "
                    "complete logistic refit"
                ),
                "go_redundancy_filter": (
                    "ancestor/descendant-related GO terms with >=0.65 Jaccard "
                    "overlap in significant terminal units are not both selected"
                ),
            },
            "inputs": [
                {
                    "basename": path.name,
                    "sha256": sha256(path),
                }
                for path in (
                    analysis_manifest,
                    model_path,
                    all_path,
                    args.go_obo,
                    args.kegg_pathway_names,
                )
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in outputs
            ],
            "software": {
                "python": os.sys.version.split()[0],
                "matplotlib": matplotlib.__version__,
            },
        }
        (temporary / "figure_manifest.json").write_text(
            json.dumps(figure_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
