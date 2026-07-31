#!/usr/bin/env python3
"""Compare paired gene-loss calls among genome units from the same species."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import binomtest, chi2


RESOLVED = {"retained", "decayed", "deleted"}
POSITIVE = {"decayed", "deleted"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_metadata(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, list[str]]]:
    rows: dict[str, dict[str, str]] = {}
    groups: dict[str, list[str]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("include", "").lower() != "true":
                continue
            unit = row["assembly_unit_id"]
            rows[unit] = row
            groups[row["biological_species"]].append(unit)
    groups = {species: units for species, units in groups.items() if len(units) >= 2}
    return rows, groups


def read_matrix(
    path: Path,
    selected_units: set[str],
) -> dict[str, dict[str, str]]:
    calls: dict[str, dict[str, str]] = {unit: {} for unit in selected_units}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "reference_gene_id",
            "assembly_unit_id",
            "manuscript_classification",
            "manuscript_positive_loss",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Matrix is missing columns: {sorted(missing)}")
        for row in reader:
            unit = row["assembly_unit_id"]
            if unit not in selected_units:
                continue
            gene = row["reference_gene_id"]
            state = row["manuscript_classification"]
            expected_positive = state in POSITIVE
            observed_positive = row["manuscript_positive_loss"].lower() == "true"
            if state in RESOLVED and observed_positive != expected_positive:
                raise ValueError(f"Inconsistent positive flag for {unit}, {gene}")
            if gene in calls[unit]:
                raise ValueError(f"Duplicate row for {unit}, {gene}")
            calls[unit][gene] = state
    return calls


def bh_adjust(pvalues: list[float]) -> list[float]:
    count = len(pvalues)
    order = sorted(range(count), key=pvalues.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        value = min(1.0, pvalues[original_index] * count / rank)
        running = min(running, value)
        adjusted[original_index] = running
    return adjusted


def exact_mcnemar(b: int, c: int) -> tuple[float, float, float, float]:
    discordant = b + c
    if discordant == 0:
        return 1.0, 1.0, 0.0, math.inf
    test = binomtest(b, discordant, 0.5, alternative="two-sided")
    interval = test.proportion_ci(confidence_level=0.95, method="exact")

    def odds(probability: float) -> float:
        if probability <= 0:
            return 0.0
        if probability >= 1:
            return math.inf
        return probability / (1.0 - probability)

    matched_or = b / c if c else math.inf
    return test.pvalue, matched_or, odds(interval.low), odds(interval.high)


def cochran_q(binary_rows: list[list[int]]) -> tuple[float, int, float]:
    if not binary_rows:
        raise ValueError("Cochran's Q requires at least one matched gene")
    k = len(binary_rows[0])
    column_sums = [sum(row[column] for row in binary_rows) for column in range(k)]
    row_sums = [sum(row) for row in binary_rows]
    total = sum(column_sums)
    denominator = k * total - sum(value * value for value in row_sums)
    if denominator == 0:
        return 0.0, k - 1, 1.0
    numerator = (k - 1) * (k * sum(value * value for value in column_sums) - total * total)
    statistic = numerator / denominator
    return statistic, k - 1, float(chi2.sf(statistic, k - 1))


def format_number(value: float | int | str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if math.isinf(value):
        return "Inf"
    return f"{value:.12g}"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_number(row.get(key, "")) for key in fieldnames})


def main() -> None:
    args = parse_args()
    metadata, groups = read_metadata(args.metadata)
    selected_units = {unit for units in groups.values() for unit in units}
    calls = read_matrix(args.matrix, selected_units)
    args.output_dir.mkdir(parents=True, exist_ok=False)

    unit_rows: list[dict[str, object]] = []
    omnibus_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []

    for species in sorted(groups):
        units = groups[species]
        common_genes = set.intersection(
            *[
                {gene for gene, state in calls[unit].items() if state in RESOLVED}
                for unit in units
            ]
        )
        ordered_genes = sorted(common_genes)
        binary_rows = [
            [int(calls[unit][gene] in POSITIVE) for unit in units]
            for gene in ordered_genes
        ]

        for unit_index, unit in enumerate(units):
            states = [calls[unit][gene] for gene in ordered_genes]
            counts = {state: states.count(state) for state in sorted(RESOLVED)}
            positive = counts["decayed"] + counts["deleted"]
            interval = binomtest(positive, len(states)).proportion_ci(
                confidence_level=0.95,
                method="wilson",
            )
            unit_rows.append(
                {
                    "biological_species": species,
                    "assembly_unit_id": unit,
                    "unit_label": metadata[unit]["haplotype_or_subgenome"],
                    "matched_resolved_genes": len(states),
                    "retained": counts["retained"],
                    "decayed": counts["decayed"],
                    "deleted": counts["deleted"],
                    "positive_decayed_plus_deleted": positive,
                    "loss_rate": positive / len(states),
                    "loss_rate_percent": 100.0 * positive / len(states),
                    "wilson_95ci_low": interval.low,
                    "wilson_95ci_high": interval.high,
                }
            )

        if len(units) > 2:
            statistic, df, pvalue = cochran_q(binary_rows)
            test_name = "Cochran_Q"
        else:
            b = sum(row[0] == 1 and row[1] == 0 for row in binary_rows)
            c = sum(row[0] == 0 and row[1] == 1 for row in binary_rows)
            pvalue, _, _, _ = exact_mcnemar(b, c)
            statistic, df = "", ""
            test_name = "exact_McNemar"
        omnibus_rows.append(
            {
                "biological_species": species,
                "genome_unit_count": len(units),
                "matched_resolved_genes": len(ordered_genes),
                "test": test_name,
                "statistic": statistic,
                "df": df,
                "p_value": pvalue,
                "significant_0.05": str(pvalue < 0.05).lower(),
            }
        )

        species_pairwise_start = len(pairwise_rows)
        for first_index in range(len(units)):
            for second_index in range(first_index + 1, len(units)):
                first = units[first_index]
                second = units[second_index]
                b = sum(
                    row[first_index] == 1 and row[second_index] == 0
                    for row in binary_rows
                )
                c = sum(
                    row[first_index] == 0 and row[second_index] == 1
                    for row in binary_rows
                )
                both_loss = sum(
                    row[first_index] == 1 and row[second_index] == 1
                    for row in binary_rows
                )
                neither_loss = len(binary_rows) - both_loss - b - c
                pvalue, matched_or, ci_low, ci_high = exact_mcnemar(b, c)
                first_rate = sum(row[first_index] for row in binary_rows) / len(binary_rows)
                second_rate = sum(row[second_index] for row in binary_rows) / len(binary_rows)
                pairwise_rows.append(
                    {
                        "biological_species": species,
                        "unit_1": first,
                        "unit_1_label": metadata[first]["haplotype_or_subgenome"],
                        "unit_2": second,
                        "unit_2_label": metadata[second]["haplotype_or_subgenome"],
                        "matched_resolved_genes": len(binary_rows),
                        "both_loss": both_loss,
                        "unit_1_loss_only": b,
                        "unit_2_loss_only": c,
                        "neither_loss": neither_loss,
                        "unit_1_loss_rate": first_rate,
                        "unit_2_loss_rate": second_rate,
                        "rate_difference_percentage_points": 100.0 * (first_rate - second_rate),
                        "matched_odds_ratio": matched_or,
                        "matched_or_95ci_low": ci_low,
                        "matched_or_95ci_high": ci_high,
                        "exact_mcnemar_p": pvalue,
                    }
                )
        species_rows = pairwise_rows[species_pairwise_start:]
        adjusted = bh_adjust([float(row["exact_mcnemar_p"]) for row in species_rows])
        for row, qvalue in zip(species_rows, adjusted):
            row["bh_adjusted_p_within_species"] = qvalue
            row["significant_bh_0.05"] = str(qvalue < 0.05).lower()

    write_tsv(
        args.output_dir / "unit_loss_rates_matched.tsv",
        [
            "biological_species",
            "assembly_unit_id",
            "unit_label",
            "matched_resolved_genes",
            "retained",
            "decayed",
            "deleted",
            "positive_decayed_plus_deleted",
            "loss_rate",
            "loss_rate_percent",
            "wilson_95ci_low",
            "wilson_95ci_high",
        ],
        unit_rows,
    )
    write_tsv(
        args.output_dir / "within_species_omnibus_tests.tsv",
        [
            "biological_species",
            "genome_unit_count",
            "matched_resolved_genes",
            "test",
            "statistic",
            "df",
            "p_value",
            "significant_0.05",
        ],
        omnibus_rows,
    )
    write_tsv(
        args.output_dir / "within_species_pairwise_mcnemar.tsv",
        [
            "biological_species",
            "unit_1",
            "unit_1_label",
            "unit_2",
            "unit_2_label",
            "matched_resolved_genes",
            "both_loss",
            "unit_1_loss_only",
            "unit_2_loss_only",
            "neither_loss",
            "unit_1_loss_rate",
            "unit_2_loss_rate",
            "rate_difference_percentage_points",
            "matched_odds_ratio",
            "matched_or_95ci_low",
            "matched_or_95ci_high",
            "exact_mcnemar_p",
            "bh_adjusted_p_within_species",
            "significant_bh_0.05",
        ],
        pairwise_rows,
    )

    manifest = {
        "status": "PASS_WITHIN_SPECIES_GENOME_UNIT_LOSS_COMPARISON",
        "positive_definition": "manuscript decayed + deleted",
        "comparison_denominator": (
            "reference genes classified as retained, decayed, or deleted in every "
            "genome unit of the focal species"
        ),
        "omnibus_test": "Cochran's Q for species with three or more genome units",
        "two_unit_test": "two-sided exact McNemar test",
        "pairwise_test": "two-sided exact McNemar test",
        "multiple_testing": "Benjamini-Hochberg within each species",
        "matrix_sha256": sha256(args.matrix),
        "metadata_sha256": sha256(args.metadata),
        "species_count": len(groups),
        "genome_unit_count": len(selected_units),
        "pairwise_comparison_count": len(pairwise_rows),
    }
    (args.output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
