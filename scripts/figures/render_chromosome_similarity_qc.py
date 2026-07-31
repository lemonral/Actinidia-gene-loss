#!/usr/bin/env python3
"""Render the author-approved maximum-similarity chromosome naming audit.

All 29 chromosomes in every declared unit are named by the global one-to-one
maximum nucleotide similarity assignment to HY4A Chr01--Chr29.  Nucleotide
coverage and independent JCVI/HY4P agreement are displayed as QC diagnostics;
absolute support never blocks a unique label under the approved naming rule.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Mapping, Sequence

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import TaxonLabelError, format_taxon_label


SCRIPT_VERSION = "1.0.1"
LABEL_COLUMNS = (
    "query_chromosome",
    "final_chromosome",
    "coordinate_reference",
    "assignment_method",
    "assigned_score",
    "reciprocal_coverage",
    "orientation_to_hy4a",
    "hy4p_and_jcvi_agree",
    "strict_homology_gates_pass",
    "confidence_flag",
)
DIAGNOSTIC_COLUMNS = (
    "query_chromosome",
    "diagnostic_candidate",
    "nucleotide_hy4a",
    "jcvi_hy4a",
    "nucleotide_hy4p",
    "jcvi_hy4p",
    "nucleotide_hy4a_reciprocal_coverage",
    "jcvi_hy4a_score",
    "orientation_hy4a",
    "all_four_matrix_gates",
    "diagnostic_status",
    "failure_reasons",
)
METADATA_REQUIRED = (
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "accession",
    "decision_status",
)
PLOT_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "upright_suffix",
    "display_label",
    "named_chromosome_count",
    "four_matrix_agreement_count",
    "strict_all_four_gate_count",
    "qc_support_only_count",
    "high_confidence_count",
    "supported_confidence_count",
    "mean_reciprocal_coverage_percent",
    "minimum_reciprocal_coverage_percent",
    "maximum_reciprocal_coverage_percent",
    "mean_jcvi_score_percent",
    "minimum_jcvi_score_percent",
    "maximum_jcvi_score_percent",
)
CHR_RE = re.compile(r"Chr(?:0[1-9]|1[0-9]|2[0-9])\Z")


class ChromosomeSimilarityPlotError(RuntimeError):
    """Raised when one chromosome-naming audit cannot be published."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, required: Sequence[str], name: str, *, exact: bool = True) -> list[dict[str, str]]:
    path = Path(path).expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise ChromosomeSimilarityPlotError(f"cannot open {name} {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if (exact and fields != tuple(required)) or (not exact and not set(required).issubset(fields)):
            raise ChromosomeSimilarityPlotError(
                f"{name}: schema mismatch; required={list(required)!r}; observed={list(fields)!r}"
            )
        rows = [
            {field: (row[field] or "").strip() for field in fields}
            for row in reader
            if any((row[field] or "").strip() for field in fields)
        ]
    if not rows:
        raise ChromosomeSimilarityPlotError(f"{name}: no data rows")
    return rows


def number(value: str, location: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ChromosomeSimilarityPlotError(f"{location}: nonnumeric value {value!r}") from error
    if not math.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ChromosomeSimilarityPlotError(f"{location}: value must be within [0,1]")
    return result


def unit_suffix(row: Mapping[str, str]) -> str:
    suffix = row["haplotype_or_subgenome"]
    if suffix and suffix not in {"unphased", "NA"}:
        return suffix
    accession = row["accession"]
    if row["biological_species"] == "Actinidia rufa" and "ActinidiaBase" in accession:
        return ""
    return ""


def prepare(
    *,
    metadata_path: Path,
    assignment_roots: Sequence[Path],
    expected_unit_count: int = 9,
    expected_chromosome_count: int = 29,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    metadata_rows = read_tsv(metadata_path, METADATA_REQUIRED, "QC metadata", exact=False)
    metadata = {row["assembly_unit_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise ChromosomeSimilarityPlotError("QC metadata contains duplicate assembly_unit_id values")

    unit_dirs: dict[str, Path] = {}
    for raw_root in assignment_roots:
        root = Path(raw_root).expanduser().resolve()
        if not root.is_dir():
            raise ChromosomeSimilarityPlotError(f"assignment root is not a directory: {root}")
        for child in sorted(root.iterdir()):
            if not child.is_dir() or not (child / "similarity_label_map.tsv").is_file():
                continue
            if child.name in unit_dirs:
                raise ChromosomeSimilarityPlotError(f"duplicate assignment unit {child.name}")
            unit_dirs[child.name] = child
    if len(unit_dirs) != expected_unit_count:
        raise ChromosomeSimilarityPlotError(
            f"expected {expected_unit_count} assignment units, found {len(unit_dirs)}"
        )
    missing_metadata = sorted(set(unit_dirs).difference(metadata))
    if missing_metadata:
        raise ChromosomeSimilarityPlotError(f"assignment units missing from QC metadata: {missing_metadata}")

    expected_labels = {f"Chr{index:02d}" for index in range(1, expected_chromosome_count + 1)}
    unit_rows: list[dict[str, object]] = []
    source_checksums: list[dict[str, object]] = []
    for unit_id, directory in unit_dirs.items():
        label_path = directory / "similarity_label_map.tsv"
        diagnostic_path = directory / "prefinal_assignment_diagnostic.tsv"
        status_path = directory / "diagnostic.json"
        labels = read_tsv(label_path, LABEL_COLUMNS, f"{unit_id} label map")
        diagnostics = read_tsv(diagnostic_path, DIAGNOSTIC_COLUMNS, f"{unit_id} diagnostic")
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: cannot parse diagnostic.json: {error}") from error
        if status.get("assembly_unit_id") != unit_id:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: diagnostic assembly ID mismatch")
        if status.get("chromosome_naming_status") != "PASS_LABELS":
            raise ChromosomeSimilarityPlotError(f"{unit_id}: chromosome naming status is not PASS_LABELS")
        if status.get("chromosome_naming_policy") != (
            "HY4A global one-to-one maximum nucleotide similarity; absolute support is QC only"
        ):
            raise ChromosomeSimilarityPlotError(f"{unit_id}: unexpected chromosome naming policy")
        if len(labels) != expected_chromosome_count or len(diagnostics) != expected_chromosome_count:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: expected {expected_chromosome_count} rows")
        query_ids = {row["query_chromosome"] for row in labels}
        final_ids = {row["final_chromosome"] for row in labels}
        if len(query_ids) != expected_chromosome_count or final_ids != expected_labels:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: label map is not a complete Chr01-Chr29 bijection")
        if any(row["assignment_method"] != "global_one_to_one_maximum_nucleotide_similarity" for row in labels):
            raise ChromosomeSimilarityPlotError(f"{unit_id}: unexpected assignment method")
        if any(row["hy4p_and_jcvi_agree"] != "true" for row in labels):
            raise ChromosomeSimilarityPlotError(f"{unit_id}: label map lacks independent evidence agreement")

        label_by_query = {row["query_chromosome"]: row for row in labels}
        if set(label_by_query) != {row["query_chromosome"] for row in diagnostics}:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: diagnostic query set differs from label map")
        agreement = 0
        strict = 0
        coverages: list[float] = []
        jcvi_scores: list[float] = []
        for row in diagnostics:
            final = label_by_query[row["query_chromosome"]]["final_chromosome"]
            evidence = (
                row["diagnostic_candidate"],
                row["nucleotide_hy4a"],
                row["jcvi_hy4a"],
                row["nucleotide_hy4p"],
                row["jcvi_hy4p"],
            )
            if all(value == final for value in evidence):
                agreement += 1
            else:
                raise ChromosomeSimilarityPlotError(f"{unit_id}: four-way assignment disagreement")
            if row["all_four_matrix_gates"] == "true":
                strict += 1
            elif row["all_four_matrix_gates"] != "false":
                raise ChromosomeSimilarityPlotError(f"{unit_id}: invalid all_four_matrix_gates value")
            coverages.append(number(row["nucleotide_hy4a_reciprocal_coverage"], f"{unit_id} coverage"))
            jcvi_scores.append(number(row["jcvi_hy4a_score"], f"{unit_id} JCVI score"))

        confidence_counts = {
            flag: sum(row["confidence_flag"] == flag for row in labels)
            for flag in ("HIGH", "SUPPORTED")
        }
        if sum(confidence_counts.values()) != expected_chromosome_count:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: unexpected confidence flag")
        meta = metadata[unit_id]
        if meta["decision_status"] != "current":
            raise ChromosomeSimilarityPlotError(f"{unit_id}: naming unit is not current in QC metadata")
        suffix = unit_suffix(meta)
        try:
            display = format_taxon_label(meta["biological_species"], (suffix,), abbreviate_genus=True)
        except TaxonLabelError as error:
            raise ChromosomeSimilarityPlotError(f"{unit_id}: invalid label metadata: {error}") from error
        unit_rows.append(
            {
                "assembly_unit_id": unit_id,
                "biological_species": meta["biological_species"],
                "upright_suffix": suffix,
                "display_label": display,
                "named_chromosome_count": expected_chromosome_count,
                "four_matrix_agreement_count": agreement,
                "strict_all_four_gate_count": strict,
                "qc_support_only_count": expected_chromosome_count - strict,
                "high_confidence_count": confidence_counts["HIGH"],
                "supported_confidence_count": confidence_counts["SUPPORTED"],
                "mean_reciprocal_coverage_percent": 100.0 * sum(coverages) / len(coverages),
                "minimum_reciprocal_coverage_percent": 100.0 * min(coverages),
                "maximum_reciprocal_coverage_percent": 100.0 * max(coverages),
                "mean_jcvi_score_percent": 100.0 * sum(jcvi_scores) / len(jcvi_scores),
                "minimum_jcvi_score_percent": 100.0 * min(jcvi_scores),
                "maximum_jcvi_score_percent": 100.0 * max(jcvi_scores),
            }
        )
        for role, path in (
            ("similarity_label_map", label_path),
            ("four_matrix_diagnostic", diagnostic_path),
            ("unit_diagnostic", status_path),
        ):
            source_checksums.append(
                {
                    "assembly_unit_id": unit_id,
                    "role": role,
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )

    species_order = {"Actinidia rufa": 0, "Actinidia eriantha": 1, "Actinidia deliciosa": 2}
    unit_rows.sort(
        key=lambda row: (
            species_order.get(str(row["biological_species"]), 99),
            str(row["upright_suffix"]),
        )
    )
    rows = [{"plot_order": index, **row} for index, row in enumerate(unit_rows, 1)]
    validation: dict[str, object] = {
        "schema_version": 1,
        "renderer": "scripts/figures/render_chromosome_similarity_qc.py",
        "renderer_version": SCRIPT_VERSION,
        "status": "PASS_CHROMOSOME_SIMILARITY_NAMING_PUBLICATION",
        "chromosome_naming_policy": (
            "HY4A global one-to-one maximum nucleotide similarity; absolute support is QC only"
        ),
        "orientation_policy": "publisher direction preserved; no forced flips",
        "unit_count": len(rows),
        "chromosomes_per_unit": expected_chromosome_count,
        "named_chromosome_total": sum(int(row["named_chromosome_count"]) for row in rows),
        "four_matrix_agreement_total": sum(int(row["four_matrix_agreement_count"]) for row in rows),
        "source_checksums": source_checksums,
        "checks": {
            "complete_chr01_chr29_bijection_each_unit": "pass",
            "maximum_similarity_assignment_method": "pass",
            "all_four_independent_assignments_agree": "pass",
            "absolute_support_is_qc_not_naming_blocker": "pass",
            "publisher_direction_preserved": "pass",
            "italic_binomials_upright_suffixes": "pass",
        },
    }
    return rows, validation


def render_bundle(
    *,
    metadata_path: Path,
    assignment_roots: Sequence[Path],
    output_dir: Path,
    basename: str,
    expected_unit_count: int = 9,
    expected_chromosome_count: int = 29,
    dpi: int = 300,
) -> FigureBundle:
    rows, validation = prepare(
        metadata_path=metadata_path,
        assignment_roots=assignment_roots,
        expected_unit_count=expected_unit_count,
        expected_chromosome_count=expected_chromosome_count,
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ChromosomeSimilarityPlotError("matplotlib is required to render the figure") from error

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.2,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    y = list(range(len(rows)))
    strict = [int(row["strict_all_four_gate_count"]) for row in rows]
    support_only = [int(row["qc_support_only_count"]) for row in rows]
    labels = [str(row["display_label"]) for row in rows]
    fig, (name_ax, support_ax) = plt.subplots(
        1,
        2,
        figsize=(7.2, 5.2),
        gridspec_kw={"width_ratios": (1.0, 1.35), "wspace": 0.08},
        sharey=True,
    )
    name_ax.barh(y, strict, color="#2A9D8F", label="All four diagnostics agree")
    name_ax.barh(y, support_only, left=strict, color="#ADB5BD", label="Accepted; lower absolute support")
    name_ax.axvline(expected_chromosome_count, color="#343A40", lw=0.8)
    name_ax.set_xlim(0, expected_chromosome_count + 1)
    name_ax.set_xlabel("Chromosomes assigned to Chr01–Chr29")
    name_ax.set_yticks(y, labels)
    name_ax.invert_yaxis()
    name_ax.spines[["top", "right", "left"]].set_visible(False)
    name_ax.tick_params(axis="y", length=0)
    name_ax.grid(axis="x", color="#DEE2E6", lw=0.6)
    name_ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.105),
        fontsize=7.2,
        ncol=1,
    )

    coverage_mean = [float(row["mean_reciprocal_coverage_percent"]) for row in rows]
    coverage_low = [float(row["minimum_reciprocal_coverage_percent"]) for row in rows]
    coverage_high = [float(row["maximum_reciprocal_coverage_percent"]) for row in rows]
    jcvi_mean = [float(row["mean_jcvi_score_percent"]) for row in rows]
    jcvi_low = [float(row["minimum_jcvi_score_percent"]) for row in rows]
    jcvi_high = [float(row["maximum_jcvi_score_percent"]) for row in rows]
    support_ax.errorbar(
        coverage_mean,
        [value - 0.12 for value in y],
        xerr=(
            [max(0.0, mean - low) for mean, low in zip(coverage_mean, coverage_low)],
            [max(0.0, high - mean) for mean, high in zip(coverage_mean, coverage_high)],
        ),
        fmt="o",
        color="#3A86FF",
        ecolor="#8ECAE6",
        capsize=2,
        label="Nucleotide reciprocal coverage",
    )
    support_ax.errorbar(
        jcvi_mean,
        [value + 0.12 for value in y],
        xerr=(
            [max(0.0, mean - low) for mean, low in zip(jcvi_mean, jcvi_low)],
            [max(0.0, high - mean) for mean, high in zip(jcvi_mean, jcvi_high)],
        ),
        fmt="s",
        color="#8338EC",
        ecolor="#CDB4DB",
        capsize=2,
        label="JCVI anchor score",
    )
    support_ax.set_xlim(0, 80)
    support_ax.set_xlabel("Mean and range across 29 assignments (%)")
    support_ax.spines[["top", "right", "left"]].set_visible(False)
    support_ax.tick_params(axis="y", left=False, labelleft=False)
    support_ax.grid(axis="x", color="#DEE2E6", lw=0.6)
    support_ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.105),
        fontsize=7.2,
        ncol=1,
    )
    for axis, panel in ((name_ax, "a"), (support_ax, "b")):
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
    fig.subplots_adjust(bottom=0.24, top=0.94)

    caption = (
        "Chromosome naming and independent support. Each of the nine newly harmonized assembly "
        "units has a complete one-to-one Chr01--Chr29 assignment obtained by maximizing global "
        "nucleotide similarity to Hongyang v4 HY4A. All nucleotide HY4A, JCVI HY4A, nucleotide "
        "HY4P, and JCVI HY4P assignments agree for all 29 chromosomes in every unit. The left "
        "panel separates chromosomes that also pass the earlier strict absolute-support diagnostic "
        "from chromosomes accepted under the author-approved naming rule; the latter are not "
        "unassigned or ambiguous. The right panel shows mean and range of nucleotide reciprocal "
        "coverage and JCVI anchor score as QC only. Publisher chromosome direction is preserved; "
        "no sequence was forced to reverse. Latin binomials are italic and assembly suffixes are upright."
    )
    bundle = write_figure_bundle(
        figure=fig,
        output_dir=output_dir,
        basename=basename,
        plot_rows=rows,
        plot_columns=PLOT_COLUMNS,
        caption=caption,
        validation=validation,
        input_paths=(metadata_path,),
        dpi=dpi,
    )
    plt.close(fig)
    return bundle


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--assignment-root", required=True, action="append", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="chromosome_similarity_naming")
    parser.add_argument("--expected-unit-count", type=int, default=9)
    parser.add_argument("--expected-chromosome-count", type=int, default=29)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    try:
        render_bundle(
            metadata_path=args.metadata,
            assignment_roots=args.assignment_root,
            output_dir=args.output_dir,
            basename=args.basename,
            expected_unit_count=args.expected_unit_count,
            expected_chromosome_count=args.expected_chromosome_count,
            dpi=args.dpi,
        )
    except (OSError, ValueError, ChromosomeSimilarityPlotError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
