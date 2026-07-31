#!/usr/bin/env python3
"""Analyze chromosome-scale positions of manuscript-method decayed loci.

Only article-method ``decayed`` calls with an observed target-assembly residual
coordinate enter spatial numerators.  Target-assembly GFF gene features provide
the opportunity denominator.  ``deleted`` calls never enter this analysis.
"""

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
from collections import Counter
from pathlib import Path


STRICT_GROUPS = (
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
)
ANALYSIS_GROUPS = (
    "all_decayed",
    "strict_pseudogenized",
    "non_strict_decayed",
    *STRICT_GROUPS,
)
ZONE_LABELS = (
    "Z1_terminal",
    "Z2_subterminal",
    "Z3_intermediate_outer",
    "Z4_intermediate_inner",
    "Z5_central",
)


class DecayedPositionError(ValueError):
    """Raised when frozen coordinate evidence cannot be reconciled."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-positions", required=True, type=Path)
    parser.add_argument("--gff-registry", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-chromosomes", type=int, default=29)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise DecayedPositionError(f"missing or empty input: {path}")
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    binding(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = list(reader.fieldnames or [])
    if not rows or not fields or len(fields) != len(set(fields)):
        raise DecayedPositionError(f"invalid TSV: {path.name}")
    return rows, fields


def require(fields: list[str], expected: set[str], label: str) -> None:
    missing = sorted(expected - set(fields))
    if missing:
        raise DecayedPositionError(
            f"{label} missing fields: {', '.join(missing)}"
        )


def resolve(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise DecayedPositionError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise DecayedPositionError(f"path escapes data root: {value!r}")
    return path


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            import io

            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                writer.writerows(rows)


def parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key] = value
    return attributes


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current = ""
    length = 0
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if current:
                    lengths[current] = length
                current = raw[1:].strip().split()[0]
                if not current or current in lengths:
                    raise DecayedPositionError(
                        f"{path.name}:{line_number}: invalid FASTA identifier"
                    )
                length = 0
            else:
                sequence = raw.strip()
                if not current or not sequence:
                    if not sequence:
                        continue
                    raise DecayedPositionError(
                        f"{path.name}:{line_number}: sequence before header"
                    )
                length += len(sequence)
    if current:
        lengths[current] = length
    if not lengths or any(value < 1 for value in lengths.values()):
        raise DecayedPositionError(f"{path.name}: invalid FASTA")
    return lengths


def zone_for_distance(value: float) -> str:
    if not 0.0 <= value <= 1.0 + 1e-12:
        raise DecayedPositionError(f"normalized end distance outside [0,1]: {value}")
    index = min(int(min(value, 1.0) * 5.0), 4)
    return ZONE_LABELS[index]


def normalized_end_distance(midpoint: float, length: int) -> float:
    if length < 2 or not 1.0 <= midpoint <= length:
        raise DecayedPositionError(
            f"invalid midpoint/length combination: {midpoint}/{length}"
        )
    return min(midpoint - 1.0, length - midpoint) / ((length - 1.0) / 2.0)


def poisson_interval(count: int, opportunity: int) -> tuple[float, float]:
    from scipy.stats import chi2

    if opportunity <= 0:
        raise DecayedPositionError("non-positive gene opportunity")
    lower_count = 0.0 if count == 0 else 0.5 * chi2.ppf(0.025, 2 * count)
    upper_count = 0.5 * chi2.ppf(0.975, 2 * (count + 1))
    return (
        lower_count / opportunity * 1000.0,
        upper_count / opportunity * 1000.0,
    )


def burden_row(count: int, genes: int) -> dict[str, object]:
    lower, upper = poisson_interval(count, genes)
    return {
        "decayed_loci": count,
        "target_annotated_genes": genes,
        "decayed_loci_per_1000_genes": f"{count / genes * 1000.0:.12g}",
        "poisson_95ci_lower_per_1000": f"{lower:.12g}",
        "poisson_95ci_upper_per_1000": f"{upper:.12g}",
    }


def bh_adjust(rows: list[dict[str, object]], key: str, output: str) -> None:
    ordered = sorted(enumerate(rows), key=lambda item: float(item[1][key]))
    count = len(ordered)
    running = 1.0
    for reverse_rank, (index, row) in enumerate(reversed(ordered), 1):
        rank = count - reverse_rank + 1
        adjusted = min(running, float(row[key]) * count / rank, 1.0)
        running = adjusted
        rows[index][output] = f"{adjusted:.12g}"


def opportunity_chi_square(
    observed: list[int],
    opportunities: list[int],
) -> tuple[float, int, float]:
    from scipy.stats import chi2

    if len(observed) != len(opportunities) or len(observed) < 2:
        raise DecayedPositionError("invalid opportunity chi-square input")
    total_observed = sum(observed)
    total_opportunity = sum(opportunities)
    if total_observed == 0 or total_opportunity == 0:
        return 0.0, len(observed) - 1, 1.0
    expected = [
        total_observed * opportunity / total_opportunity
        for opportunity in opportunities
    ]
    statistic = sum(
        (actual - model) ** 2 / model
        for actual, model in zip(observed, expected, strict=True)
        if model > 0
    )
    degrees = len(observed) - 1
    return statistic, degrees, float(chi2.sf(statistic, degrees))


def fit_negative_binomial(
    rows: list[dict[str, object]],
    *,
    effect: str,
    levels: list[str],
    adjustment_terms: list[str],
    reference_level: str | None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    import numpy as np
    import pandas as pd
    import patsy
    import statsmodels.formula.api as smf
    from scipy.stats import chi2, norm

    data = pd.DataFrame(rows)
    if (data["target_annotated_genes"].astype(int) <= 0).any():
        raise DecayedPositionError("model contains a non-positive denominator")
    offset = np.log(data["target_annotated_genes"].astype(float))
    adjustment = " + ".join(f"C({term})" for term in adjustment_terms)
    effect_term = f"C({effect})"
    full_formula = f"decayed_loci ~ {effect_term} + {adjustment}"
    reduced_formula = f"decayed_loci ~ {adjustment}"
    full = smf.negativebinomial(
        full_formula,
        data=data,
        offset=offset,
    ).fit(method="bfgs", maxiter=1000, disp=False)
    reduced = smf.negativebinomial(
        reduced_formula,
        data=data,
        offset=offset,
    ).fit(method="bfgs", maxiter=1000, disp=False)
    if not bool(full.mle_retvals.get("converged")) or not bool(
        reduced.mle_retvals.get("converged")
    ):
        raise DecayedPositionError(f"negative-binomial {effect} model did not converge")

    statistic = max(0.0, 2.0 * (float(full.llf) - float(reduced.llf)))
    degrees = int(len(full.params) - len(reduced.params))
    p_value = float(chi2.sf(statistic, degrees))
    alpha = float(full.params["alpha"])
    summary = {
        "effect": effect,
        "model": "negative_binomial_log_link_offset_log_target_gene_count",
        "rows": len(data),
        "total_decayed_loci": int(data["decayed_loci"].sum()),
        "total_target_annotated_genes": int(
            data["target_annotated_genes"].sum()
        ),
        "dispersion_alpha": f"{alpha:.12g}",
        "full_log_likelihood": f"{float(full.llf):.12g}",
        "reduced_log_likelihood": f"{float(reduced.llf):.12g}",
        "likelihood_ratio_chi_square": f"{statistic:.12g}",
        "degrees_of_freedom": degrees,
        "p_value": f"{p_value:.12g}",
    }

    design_info = full.model.data.design_info
    # The discrete NegativeBinomial model appends the estimated ``alpha`` to
    # ``model.exog_names`` even though it is not a Patsy design column.
    exog_names = list(design_info.column_names)
    beta = full.params.loc[exog_names].to_numpy(dtype=float)
    covariance = full.cov_params().loc[exog_names, exog_names].to_numpy(dtype=float)
    grid_columns = [effect, *adjustment_terms]
    grid_rows: list[dict[str, str]] = []
    unique_adjustments = [
        sorted(str(value) for value in data[term].unique())
        for term in adjustment_terms
    ]
    import itertools

    adjustment_grid = list(itertools.product(*unique_adjustments))
    for level in levels:
        for values in adjustment_grid:
            row = {effect: level}
            row.update(dict(zip(adjustment_terms, values, strict=True)))
            grid_rows.append(row)
    grid = pd.DataFrame(grid_rows, columns=grid_columns)
    design = np.asarray(
        patsy.build_design_matrices([design_info], grid)[0],
        dtype=float,
    )
    level_design: dict[str, np.ndarray] = {}
    block = len(adjustment_grid)
    for index, level in enumerate(levels):
        level_design[level] = design[index * block : (index + 1) * block].mean(
            axis=0
        )
    grand = np.vstack([level_design[level] for level in levels]).mean(axis=0)
    reference = level_design[reference_level] if reference_level else grand

    contrasts: list[dict[str, object]] = []
    for level in levels:
        vector = level_design[level]
        log_rate = float(vector @ beta)
        rate_se = math.sqrt(max(float(vector @ covariance @ vector), 0.0))
        contrast = vector - reference
        log_ratio = float(contrast @ beta)
        ratio_se = math.sqrt(max(float(contrast @ covariance @ contrast), 0.0))
        z_score = 0.0 if ratio_se == 0 else log_ratio / ratio_se
        contrast_p = 1.0 if ratio_se == 0 else 2.0 * float(norm.sf(abs(z_score)))
        contrasts.append(
            {
                effect: level,
                "comparison": (
                    f"versus_{reference_level}"
                    if reference_level
                    else "versus_adjusted_grand_mean"
                ),
                "adjusted_decayed_loci_per_1000_genes": f"{math.exp(log_rate) * 1000.0:.12g}",
                "adjusted_rate_95ci_lower_per_1000": f"{math.exp(log_rate - 1.96 * rate_se) * 1000.0:.12g}",
                "adjusted_rate_95ci_upper_per_1000": f"{math.exp(log_rate + 1.96 * rate_se) * 1000.0:.12g}",
                "rate_ratio": f"{math.exp(log_ratio):.12g}",
                "rate_ratio_95ci_lower": f"{math.exp(log_ratio - 1.96 * ratio_se):.12g}",
                "rate_ratio_95ci_upper": f"{math.exp(log_ratio + 1.96 * ratio_se):.12g}",
                "wald_z": f"{z_score:.12g}",
                "p_value": f"{contrast_p:.12g}",
            }
        )
    bh_adjust(contrasts, "p_value", "bh_q_value")
    return summary, contrasts


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise DecayedPositionError(
            f"refusing to overwrite output directory: {args.output_dir}"
        )
    root = args.data_root.resolve()
    registry_rows, registry_fields = read_tsv(args.gff_registry)
    require(
        registry_fields,
        {
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "source_group",
            "target_gff",
            "target_genome",
        },
        args.gff_registry.name,
    )
    registry = {row["assembly_unit_id"]: row for row in registry_rows}
    if len(registry) != args.expected_units or len(registry) != len(registry_rows):
        raise DecayedPositionError("GFF registry unit count changed")

    expected_chromosomes = [
        f"Chr{index:02d}"
        for index in range(1, args.expected_chromosomes + 1)
    ]
    lengths: dict[tuple[str, str], int] = {}
    genome_inputs: list[Path] = []
    for unit in sorted(registry):
        genome_path = resolve(root, registry[unit]["target_genome"])
        binding(genome_path)
        genome_inputs.append(genome_path)
        unit_lengths = read_fasta_lengths(genome_path)
        if set(unit_lengths) != set(expected_chromosomes):
            raise DecayedPositionError(
                f"{unit}: target genome is not the expected chromosome set"
            )
        for chromosome, length in unit_lengths.items():
            lengths[(unit, chromosome)] = length

    gene_counts: Counter[tuple[str, str, str]] = Counter()
    gff_inputs: list[Path] = []
    total_gene_ids: dict[str, set[str]] = {}
    for unit in sorted(registry):
        gff_path = resolve(root, registry[unit]["target_gff"])
        binding(gff_path)
        gff_inputs.append(gff_path)
        ids: set[str] = set()
        chromosomes: set[str] = set()
        with open_text(gff_path) as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise DecayedPositionError(
                        f"{gff_path.name}:{line_number}: not nine GFF fields"
                    )
                if fields[2] != "gene":
                    continue
                chromosome = fields[0]
                if chromosome not in expected_chromosomes:
                    raise DecayedPositionError(
                        f"{unit}: non-primary gene chromosome {chromosome!r}"
                    )
                start = int(fields[3])
                end = int(fields[4])
                if not 1 <= start <= end <= lengths[(unit, chromosome)]:
                    raise DecayedPositionError(
                        f"{unit}:{chromosome}: invalid gene coordinates"
                    )
                attributes = parse_attributes(fields[8])
                gene_id = attributes.get("ID") or attributes.get("gene_id")
                if not gene_id or gene_id in ids:
                    raise DecayedPositionError(
                        f"{unit}:{gff_path.name}:{line_number}: invalid gene ID"
                    )
                ids.add(gene_id)
                chromosomes.add(chromosome)
                midpoint = (start + end) / 2.0
                zone = zone_for_distance(
                    normalized_end_distance(
                        midpoint,
                        lengths[(unit, chromosome)],
                    )
                )
                gene_counts[(unit, chromosome, zone)] += 1
        if chromosomes != set(expected_chromosomes):
            raise DecayedPositionError(f"{unit}: GFF chromosome set changed")
        total_gene_ids[unit] = ids

    residual_rows, residual_fields = read_tsv(args.residual_positions)
    require(
        residual_fields,
        {
            "assembly_unit_id",
            "reference_gene_id",
            "primary_classification",
            "loss_type_group",
            "residual_chromosome_hy4a",
            "residual_midpoint_1based",
            "target_chromosome_length",
            "normalized_end_distance",
            "spatial_eligible",
        },
        args.residual_positions.name,
    )
    all_decayed: Counter[str] = Counter()
    placed_decayed: Counter[str] = Counter()
    strict_decayed: Counter[str] = Counter()
    placed_strict: Counter[str] = Counter()
    decayed_counts: Counter[tuple[str, str, str, str]] = Counter()
    detail_rows: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for row in residual_rows:
        if row["primary_classification"] != "decayed":
            continue
        unit = row["assembly_unit_id"]
        if unit not in registry:
            raise DecayedPositionError(f"unknown decayed unit: {unit}")
        key = (unit, row["reference_gene_id"])
        if key in seen:
            raise DecayedPositionError(f"duplicate decayed unit-gene row: {key}")
        seen.add(key)
        all_decayed[unit] += 1
        strict = row["loss_type_group"] in STRICT_GROUPS
        if strict:
            strict_decayed[unit] += 1
        if row["spatial_eligible"] != "true":
            continue
        chromosome = row["residual_chromosome_hy4a"]
        if chromosome not in expected_chromosomes:
            raise DecayedPositionError(f"{unit}: invalid residual chromosome")
        midpoint = float(row["residual_midpoint_1based"])
        length = lengths[(unit, chromosome)]
        if int(row["target_chromosome_length"]) != length:
            raise DecayedPositionError(
                f"{unit}:{chromosome}: residual/GFF length binding changed"
            )
        distance = normalized_end_distance(midpoint, length)
        if abs(distance - float(row["normalized_end_distance"])) > 1e-9:
            raise DecayedPositionError(
                f"{unit}:{row['reference_gene_id']}: distance mismatch"
            )
        zone = zone_for_distance(distance)
        groups = ["all_decayed"]
        if strict:
            groups.extend(["strict_pseudogenized", row["loss_type_group"]])
            placed_strict[unit] += 1
        else:
            groups.append("non_strict_decayed")
        for group in groups:
            decayed_counts[(group, unit, chromosome, zone)] += 1
        placed_decayed[unit] += 1
        detail_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": registry[unit]["biological_species"],
                "haplotype_or_subgenome": registry[unit][
                    "haplotype_or_subgenome"
                ],
                "source_group": registry[unit]["source_group"],
                "reference_gene_id": row["reference_gene_id"],
                "decayed_evidence_group": (
                    row["loss_type_group"]
                    if strict
                    else "non_strict_decayed"
                ),
                "strict_pseudogenized": str(strict).lower(),
                "residual_chromosome_hy4a": chromosome,
                "residual_midpoint_1based": f"{midpoint:.12g}",
                "chromosome_length_bp": length,
                "normalized_end_distance_0_end_1_center": f"{distance:.12g}",
                "five_zone": zone,
                "location_relation": row.get("location_relation", ""),
                "position_source": row.get("position_source", ""),
            }
        )

    localization_rows: list[dict[str, object]] = []
    for unit in sorted(registry):
        total = all_decayed[unit]
        placed = placed_decayed[unit]
        strict_total = strict_decayed[unit]
        strict_placed = placed_strict[unit]
        localization_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": registry[unit]["biological_species"],
                "haplotype_or_subgenome": registry[unit][
                    "haplotype_or_subgenome"
                ],
                "source_group": registry[unit]["source_group"],
                "article_method_decayed_rows": total,
                "spatially_placed_decayed_rows": placed,
                "unlocalized_decayed_rows": total - placed,
                "decayed_placement_fraction": f"{placed / total:.12g}",
                "strict_pseudogenized_decayed_rows": strict_total,
                "spatially_placed_strict_rows": strict_placed,
            }
        )

    unit_chromosome_rows: list[dict[str, object]] = []
    unit_zone_rows: list[dict[str, object]] = []
    for group in ANALYSIS_GROUPS:
        for unit in sorted(registry):
            for chromosome in expected_chromosomes:
                chromosome_count = 0
                chromosome_genes = 0
                for zone in ZONE_LABELS:
                    count = decayed_counts[(group, unit, chromosome, zone)]
                    genes = gene_counts[(unit, chromosome, zone)]
                    chromosome_count += count
                    chromosome_genes += genes
                    unit_zone_rows.append(
                        {
                            "analysis_group": group,
                            "assembly_unit_id": unit,
                            "biological_species": registry[unit][
                                "biological_species"
                            ],
                            "haplotype_or_subgenome": registry[unit][
                                "haplotype_or_subgenome"
                            ],
                            "source_group": registry[unit]["source_group"],
                            "chromosome_hy4a": chromosome,
                            "five_zone": zone,
                            **burden_row(count, genes),
                        }
                    )
                unit_chromosome_rows.append(
                    {
                        "analysis_group": group,
                        "assembly_unit_id": unit,
                        "biological_species": registry[unit][
                            "biological_species"
                        ],
                        "haplotype_or_subgenome": registry[unit][
                            "haplotype_or_subgenome"
                        ],
                        "source_group": registry[unit]["source_group"],
                        "chromosome_hy4a": chromosome,
                        **burden_row(chromosome_count, chromosome_genes),
                    }
                )

    pooled_chromosome_rows: list[dict[str, object]] = []
    pooled_zone_rows: list[dict[str, object]] = []
    pooled_chromosome_zone_rows: list[dict[str, object]] = []
    for group in ANALYSIS_GROUPS:
        for chromosome in expected_chromosomes:
            count = sum(
                decayed_counts[(group, unit, chromosome, zone)]
                for unit in registry
                for zone in ZONE_LABELS
            )
            genes = sum(
                gene_counts[(unit, chromosome, zone)]
                for unit in registry
                for zone in ZONE_LABELS
            )
            pooled_chromosome_rows.append(
                {
                    "analysis_group": group,
                    "chromosome_hy4a": chromosome,
                    **burden_row(count, genes),
                }
            )
            for zone in ZONE_LABELS:
                zone_count = sum(
                    decayed_counts[(group, unit, chromosome, zone)]
                    for unit in registry
                )
                zone_genes = sum(
                    gene_counts[(unit, chromosome, zone)]
                    for unit in registry
                )
                pooled_chromosome_zone_rows.append(
                    {
                        "analysis_group": group,
                        "chromosome_hy4a": chromosome,
                        "five_zone": zone,
                        **burden_row(zone_count, zone_genes),
                    }
                )
        for zone in ZONE_LABELS:
            count = sum(
                decayed_counts[(group, unit, chromosome, zone)]
                for unit in registry
                for chromosome in expected_chromosomes
            )
            genes = sum(
                gene_counts[(unit, chromosome, zone)]
                for unit in registry
                for chromosome in expected_chromosomes
            )
            pooled_zone_rows.append(
                {
                    "analysis_group": group,
                    "five_zone": zone,
                    **burden_row(count, genes),
                }
            )

    chromosome_tests: list[dict[str, object]] = []
    zone_tests: list[dict[str, object]] = []
    for group in ANALYSIS_GROUPS:
        chrom_group = [
            row
            for row in pooled_chromosome_rows
            if row["analysis_group"] == group
        ]
        statistic, degrees, p_value = opportunity_chi_square(
            [int(row["decayed_loci"]) for row in chrom_group],
            [int(row["target_annotated_genes"]) for row in chrom_group],
        )
        chromosome_tests.append(
            {
                "analysis_group": group,
                "test_scope": "all_29_chromosomes",
                "chi_square": f"{statistic:.12g}",
                "degrees_of_freedom": degrees,
                "p_value": f"{p_value:.12g}",
                "test_definition": (
                    "decayed counts versus target-annotated-gene opportunities"
                ),
            }
        )
        all_zone_group = [
            row for row in pooled_zone_rows if row["analysis_group"] == group
        ]
        statistic, degrees, p_value = opportunity_chi_square(
            [int(row["decayed_loci"]) for row in all_zone_group],
            [int(row["target_annotated_genes"]) for row in all_zone_group],
        )
        zone_tests.append(
            {
                "analysis_group": group,
                "test_scope": "all_chromosomes_pooled",
                "chromosome_hy4a": "all",
                "chi_square": f"{statistic:.12g}",
                "degrees_of_freedom": degrees,
                "p_value": f"{p_value:.12g}",
                "bh_q_value_across_29_chromosomes": "",
                "test_definition": (
                    "five-zone decayed counts versus target-gene opportunities"
                ),
            }
        )
        chromosome_specific: list[dict[str, object]] = []
        for chromosome in expected_chromosomes:
            selected = [
                row
                for row in pooled_chromosome_zone_rows
                if row["analysis_group"] == group
                and row["chromosome_hy4a"] == chromosome
            ]
            statistic, degrees, p_value = opportunity_chi_square(
                [int(row["decayed_loci"]) for row in selected],
                [int(row["target_annotated_genes"]) for row in selected],
            )
            chromosome_specific.append(
                {
                    "analysis_group": group,
                    "test_scope": "single_chromosome",
                    "chromosome_hy4a": chromosome,
                    "chi_square": f"{statistic:.12g}",
                    "degrees_of_freedom": degrees,
                    "p_value": f"{p_value:.12g}",
                    "test_definition": (
                        "five-zone decayed counts versus target-gene opportunities"
                    ),
                }
            )
        bh_adjust(
            chromosome_specific,
            "p_value",
            "bh_q_value_across_29_chromosomes",
        )
        zone_tests.extend(chromosome_specific)

    chromosome_model_summaries: list[dict[str, object]] = []
    chromosome_contrasts: list[dict[str, object]] = []
    zone_model_summaries: list[dict[str, object]] = []
    zone_contrasts: list[dict[str, object]] = []
    for group in ANALYSIS_GROUPS:
        chromosome_data = [
            {
                "assembly_unit_id": row["assembly_unit_id"],
                "chromosome_hy4a": row["chromosome_hy4a"],
                "decayed_loci": int(row["decayed_loci"]),
                "target_annotated_genes": int(row["target_annotated_genes"]),
            }
            for row in unit_chromosome_rows
            if row["analysis_group"] == group
        ]
        summary, contrasts = fit_negative_binomial(
            chromosome_data,
            effect="chromosome_hy4a",
            levels=expected_chromosomes,
            adjustment_terms=["assembly_unit_id"],
            reference_level=None,
        )
        summary["analysis_group"] = group
        for row in contrasts:
            row["analysis_group"] = group
        chromosome_model_summaries.append(summary)
        chromosome_contrasts.extend(contrasts)

        zone_data = [
            {
                "assembly_unit_id": row["assembly_unit_id"],
                "chromosome_hy4a": row["chromosome_hy4a"],
                "five_zone": row["five_zone"],
                "decayed_loci": int(row["decayed_loci"]),
                "target_annotated_genes": int(row["target_annotated_genes"]),
            }
            for row in unit_zone_rows
            if row["analysis_group"] == group
        ]
        summary, contrasts = fit_negative_binomial(
            zone_data,
            effect="five_zone",
            levels=list(ZONE_LABELS),
            adjustment_terms=["assembly_unit_id", "chromosome_hy4a"],
            reference_level="Z5_central",
        )
        summary["analysis_group"] = group
        for row in contrasts:
            row["analysis_group"] = group
        zone_model_summaries.append(summary)
        zone_contrasts.extend(contrasts)

    output_parent = args.output_dir.parent
    output_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.",
            dir=output_parent,
        )
    )
    try:
        outputs: list[Path] = []

        def emit(name: str, rows: list[dict[str, object]], *, gzip_output=False):
            if not rows:
                raise DecayedPositionError(f"refusing to write empty output {name}")
            path = staging / name
            if gzip_output:
                write_tsv_gz(path, list(rows[0]), rows)
            else:
                write_tsv(path, list(rows[0]), rows)
            outputs.append(path)

        emit("placed_decayed_loci.tsv.gz", detail_rows, gzip_output=True)
        emit("unit_decayed_localization_summary.tsv", localization_rows)
        emit("unit_chromosome_decayed_burden.tsv", unit_chromosome_rows)
        emit("unit_chromosome_five_zone_decayed_burden.tsv", unit_zone_rows)
        emit("pooled_chromosome_decayed_burden.tsv", pooled_chromosome_rows)
        emit("pooled_five_zone_decayed_burden.tsv", pooled_zone_rows)
        emit(
            "pooled_chromosome_five_zone_decayed_burden.tsv",
            pooled_chromosome_zone_rows,
        )
        emit("chromosome_opportunity_tests.tsv", chromosome_tests)
        emit("five_zone_opportunity_tests.tsv", zone_tests)
        emit("chromosome_model_summary.tsv", chromosome_model_summaries)
        emit("chromosome_adjusted_rate_ratios.tsv", chromosome_contrasts)
        emit("five_zone_model_summary.tsv", zone_model_summaries)
        emit("five_zone_adjusted_rate_ratios.tsv", zone_contrasts)

        manifest = {
            "status": "PASS_DECAYED_CHROMOSOME_POSITION_ANALYSIS",
            "schema_version": "1.0",
            "analysis_units": len(registry),
            "chromosomes_per_unit": args.expected_chromosomes,
            "five_zone_definition": {
                "coordinate": "normalized distance from nearest chromosome end",
                "range": "0 at either end; 1 at chromosome centre",
                "bins": list(ZONE_LABELS),
            },
            "primary_numerator": (
                "article-method decayed loci with observed target residual coordinate"
            ),
            "excluded_primary_classifications": [
                "deleted",
                "retained",
                "not_called_loss",
            ],
            "opportunity_denominator": (
                "target-assembly GFF gene features in the matching chromosome/zone"
            ),
            "burden_unit": "decayed loci per 1000 target annotated genes",
            "article_method_decayed_rows": sum(all_decayed.values()),
            "spatially_placed_decayed_rows": len(detail_rows),
            "unlocalized_decayed_rows": sum(all_decayed.values()) - len(detail_rows),
            "strict_pseudogenized_placed_rows": sum(placed_strict.values()),
            "target_annotated_gene_rows": sum(len(ids) for ids in total_gene_ids.values()),
            "analysis_groups": list(ANALYSIS_GROUPS),
            "tests": {
                "between_chromosome": (
                    "negative-binomial log-link model with assembly-unit fixed "
                    "effects and log target-gene-count offset; omnibus likelihood "
                    "ratio plus chromosome contrasts"
                ),
                "within_chromosome": (
                    "negative-binomial five-zone model adjusted for assembly unit "
                    "and chromosome, plus opportunity chi-square tests and BH "
                    "correction across 29 chromosome-specific tests"
                ),
            },
            "inputs": [
                binding(args.residual_positions),
                binding(args.gff_registry),
                *[binding(path) for path in genome_inputs],
                *[binding(path) for path in gff_inputs],
            ],
        }
        manifest_path = staging / "run_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        outputs.append(manifest_path)
        checksum_rows = [
            {
                "basename": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in outputs
        ]
        checksum_path = staging / "checksums.tsv"
        write_tsv(
            checksum_path,
            ["basename", "bytes", "sha256"],
            checksum_rows,
        )
        os.replace(staging, args.output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    run(args)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
