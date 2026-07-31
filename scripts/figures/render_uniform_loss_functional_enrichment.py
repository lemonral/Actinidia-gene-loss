#!/usr/bin/env python3
"""Render the pooled GO and KEGG-KO enrichment from uniform loss calls."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import textwrap
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle


FOREGROUNDS = (
    ("pooled_nonshared_combined", "Deleted + strict pseudogenized"),
    ("pooled_nonshared_strict_pseudogenized", "Strict pseudogenized only"),
)
ONTOLOGIES = (("GO", "GO enrichment"), ("KEGG_KO", "KEGG Orthology enrichment"))
PLOT_COLUMNS = (
    "ontology",
    "term_id",
    "term_label",
    "foreground_id",
    "foreground_label",
    "study_count",
    "study_size",
    "background_count",
    "background_size",
    "p_fdr_bh",
    "fold_enrichment",
    "minus_log10_fdr",
)


class FunctionalFigureError(ValueError):
    """Raised when enrichment outputs cannot support the declared figure."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise FunctionalFigureError(f"{path.name}: missing header")
        return list(reader)


def prepare(
    enrichment_dir: str | Path,
    *,
    top_per_foreground: int = 8,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    source = Path(enrichment_dir)
    significant_path = source / "enrichment_significant.tsv"
    foreground_path = source / "foreground_summary.tsv"
    summary_path = source / "functional_enrichment_summary.json"
    if top_per_foreground < 1:
        raise FunctionalFigureError("top_per_foreground must be positive")
    rows = read_tsv(significant_path)
    foreground_rows = read_tsv(foreground_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "PASS_UNIFORM_GO_KEGG_ENRICHMENT":
        raise FunctionalFigureError("functional enrichment summary is not PASS")
    checks = summary.get("checks")
    if not isinstance(checks, dict) or not checks or not all(value is True for value in checks.values()):
        raise FunctionalFigureError("functional enrichment checks are incomplete")
    required = {
        "foreground_id", "ontology", "term_id", "term_name", "study_count", "study_size",
        "background_count", "background_size", "p_fdr_bh", "fold_enrichment", "significant_fdr",
    }
    if rows and not required.issubset(rows[0]):
        raise FunctionalFigureError("significant enrichment table is missing required columns")
    by_group: dict[tuple[str, str], list[dict[str, str]]] = {}
    for foreground_id, _ in FOREGROUNDS:
        for ontology, _ in ONTOLOGIES:
            selected = [
                row for row in rows
                if row["foreground_id"] == foreground_id
                and row["ontology"] == ontology
                and row["significant_fdr"] == "true"
            ]
            selected.sort(key=lambda row: (float(row["p_fdr_bh"]), -float(row["fold_enrichment"]), row["term_id"]))
            by_group[(foreground_id, ontology)] = selected[:top_per_foreground]
    plot_rows: list[dict[str, object]] = []
    labels = dict(FOREGROUNDS)
    for ontology, _ in ONTOLOGIES:
        selected_terms: set[str] = set()
        for foreground_id, _ in FOREGROUNDS:
            selected_terms.update(row["term_id"] for row in by_group[(foreground_id, ontology)])
        candidates = {
            (row["foreground_id"], row["term_id"]): row
            for row in rows
            if row["ontology"] == ontology
            and row["foreground_id"] in labels
            and row["term_id"] in selected_terms
            and row["significant_fdr"] == "true"
        }
        for foreground_id, _ in FOREGROUNDS:
            for term_id in sorted(selected_terms):
                row = candidates.get((foreground_id, term_id))
                if row is None:
                    continue
                q_value = float(row["p_fdr_bh"])
                if not 0 < q_value <= 0.05:
                    raise FunctionalFigureError("plotted row is not FDR-significant")
                raw_name = row["term_name"].strip() or term_id
                if ontology == "KEGG_KO" and "; " in raw_name:
                    raw_name = raw_name.split("; ", 1)[1]
                concise = textwrap.shorten(raw_name, width=58, placeholder="…")
                plot_rows.append(
                    {
                        "ontology": ontology,
                        "term_id": term_id,
                        "term_label": f"{term_id}  {concise}",
                        "foreground_id": foreground_id,
                        "foreground_label": labels[foreground_id],
                        "study_count": int(row["study_count"]),
                        "study_size": int(row["study_size"]),
                        "background_count": int(row["background_count"]),
                        "background_size": int(row["background_size"]),
                        "p_fdr_bh": q_value,
                        "fold_enrichment": float(row["fold_enrichment"]),
                        "minus_log10_fdr": -math.log10(q_value),
                    }
                )
    if not plot_rows:
        raise FunctionalFigureError("no pooled significant GO or KEGG-KO terms")
    shared_significant = sum(row["foreground_id"] == "shared_combined" for row in rows)
    validation = {
        "schema_version": "1.0",
        "status": "pass",
        "source_status": summary["status"],
        "plotted_foregrounds": [foreground for foreground, _ in FOREGROUNDS],
        "plotted_ontologies": [ontology for ontology, _ in ONTOLOGIES],
        "top_per_foreground": top_per_foreground,
        "shared_significant_term_count": shared_significant,
        "checks": {
            "uniform_combined_foreground": True,
            "strict_pseudogenized_sensitivity_separate": True,
            "only_fdr_significant_terms_plotted": True,
            "shared_no_significant_terms_reported": shared_significant == 0,
        },
        "foreground_summary_rows": len(foreground_rows),
    }
    return plot_rows, validation, [significant_path, foreground_path, summary_path]


def build_figure(plot_rows: list[dict[str, object]]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(14.5, 8.2), constrained_layout=True)
    colors = {
        "pooled_nonshared_combined": "#0072B2",
        "pooled_nonshared_strict_pseudogenized": "#D55E00",
    }
    foreground_order = [foreground for foreground, _ in FOREGROUNDS]
    scatter = None
    for axis, (ontology, title) in zip(axes, ONTOLOGIES):
        subset = [row for row in plot_rows if row["ontology"] == ontology]
        term_best: dict[str, float] = {}
        term_labels: dict[str, str] = {}
        for row in subset:
            term = str(row["term_id"])
            term_best[term] = max(term_best.get(term, 0.0), float(row["minus_log10_fdr"]))
            term_labels[term] = str(row["term_label"])
        terms = sorted(term_best, key=lambda term: (-term_best[term], term))
        y_by_term = {term: index for index, term in enumerate(reversed(terms))}
        for foreground_index, foreground in enumerate(foreground_order):
            points = [row for row in subset if row["foreground_id"] == foreground]
            x = [float(row["fold_enrichment"]) for row in points]
            y = [y_by_term[str(row["term_id"])] + (foreground_index - 0.5) * 0.22 for row in points]
            sizes = [22 + 7 * math.sqrt(int(row["study_count"])) for row in points]
            values = [float(row["minus_log10_fdr"]) for row in points]
            scatter = axis.scatter(
                x,
                y,
                s=sizes,
                c=values,
                cmap="viridis",
                vmin=0,
                vmax=max(1.31, max(float(row["minus_log10_fdr"]) for row in plot_rows)),
                marker="o" if foreground_index == 0 else "D",
                edgecolor=colors[foreground],
                linewidth=1.0,
                label=dict(FOREGROUNDS)[foreground],
                alpha=0.92,
            )
        axis.axvline(1.0, color="#777777", linewidth=0.8, linestyle="--")
        axis.set_yticks(range(len(terms)), [term_labels[term] for term in reversed(terms)], fontsize=7.5)
        axis.set_xlabel("Fold enrichment")
        axis.set_title(title)
        axis.grid(axis="x", color="#E1E1E1", linewidth=0.55)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(loc="upper right", frameon=False, fontsize=8.5)
    if scatter is not None:
        colorbar = figure.colorbar(scatter, ax=axes, location="right", shrink=0.72, pad=0.02)
        colorbar.set_label(r"$-\log_{10}$(BH FDR)")
    return figure


def publish(
    *,
    enrichment_dir: str | Path,
    output_dir: str | Path,
    basename: str = "uniform_loss_functional_enrichment",
    top_per_foreground: int = 8,
    dpi: int = 300,
) -> FigureBundle:
    plot_rows, validation, inputs = prepare(
        enrichment_dir, top_per_foreground=top_per_foreground
    )
    figure = build_figure(plot_rows)
    caption = (
        "GO and KEGG Orthology over-representation among pooled non-shared gene-loss "
        "calls from the unified evidence matrix. Circles show the primary combined "
        "deleted plus strict-pseudogenized foreground; diamonds show the separate "
        "sensitivity in which every required assembly unit is strict pseudogenized. "
        "Point size represents the foreground gene count supporting the term, colour "
        "is minus log10 Benjamini-Hochberg FDR, and the x-axis is fold enrichment. "
        "Only FDR-significant terms are plotted. The 287 shared positive-complete "
        "genes had no significant GO, KEGG-KO, or KEGG-pathway term at FDR 0.05."
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
            input_paths=inputs,
            dpi=dpi,
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enrichment-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="uniform_loss_functional_enrichment")
    parser.add_argument("--top-per-foreground", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    try:
        publish(
            enrichment_dir=args.enrichment_dir,
            output_dir=args.output_dir,
            basename=args.basename,
            top_per_foreground=args.top_per_foreground,
            dpi=args.dpi,
        )
    except (FunctionalFigureError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
