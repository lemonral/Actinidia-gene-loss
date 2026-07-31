"""Explicit statistical analyses replacing date-stamped legacy scripts.

The functions intentionally separate biological species from individual
haplotypes/subgenomes.  Haplotype-level tests are allowed as a sensitivity
analysis, but species-level aggregation is the default to avoid pseudoreplication.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from pathlib import Path

from .io_utils import SchemaError, bh_adjust, format_number, read_tsv, write_tsv


def _scipy_stats():
    try:
        from scipy import stats
    except ImportError as exc:
        raise RuntimeError(
            "This statistical command needs SciPy. Install optional dependencies with: "
            "python -m pip install -e '.[statistics]'"
        ) from exc
    return stats


def _metadata_by_sample(path: str | Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path, required=["sample_id", "biological_species", "ploidy"])
    mapping: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row["sample_id"]
        if sample in mapping:
            raise SchemaError(f"{path}: duplicate sample_id {sample!r}")
        mapping[sample] = row
    return mapping


def ploidy_comparison(
    summary_path: str | Path,
    metadata_path: str | Path,
    output_dir: str | Path,
    metric: str = "assessed_loss_rate",
    unit: str = "species",
    aggregation: str = "mean",
    polyploid_labels: tuple[str, ...] = ("tetraploid", "hexaploid"),
    diploid_label: str = "diploid",
) -> dict[str, Path]:
    """Compare polyploid and diploid loss rates with a documented analysis unit."""
    if unit not in {"species", "haplotype"}:
        raise SchemaError("unit must be species or haplotype")
    if aggregation not in {"mean", "weighted"}:
        raise SchemaError("aggregation must be mean or weighted")
    stats = _scipy_stats()
    summaries = read_tsv(summary_path, required=["sample_id", metric])
    metadata = _metadata_by_sample(metadata_path)
    merged: list[dict[str, object]] = []
    for row in summaries:
        sample = row["sample_id"]
        if sample not in metadata:
            raise SchemaError(f"{summary_path}: sample_id {sample!r} absent from {metadata_path}")
        if not row[metric]:
            continue
        merged.append({**row, **metadata[sample], "metric_value": float(row[metric])})
    if not merged:
        raise SchemaError("no nonempty metric values after summary/metadata merge")
    units: list[dict[str, object]] = []
    if unit == "haplotype":
        for row in merged:
            units.append({
                "analysis_unit": row["sample_id"], "biological_species": row["biological_species"],
                "ploidy": row["ploidy"], "metric": metric, "metric_value": row["metric_value"],
                "n_haplotypes_or_subgenomes": 1, "aggregation": "none_haplotype_sensitivity",
            })
    else:
        grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
        for row in merged:
            grouped[(str(row["biological_species"]), str(row["ploidy"]))].append(row)
        for (species, ploidy), group in sorted(grouped.items()):
            if aggregation == "weighted":
                numerator_col, denominator_col = "assessed_loss_count", "reference_gene_count"
                if not all(numerator_col in row and denominator_col in row and row[numerator_col] and row[denominator_col] for row in group):
                    raise SchemaError(
                        f"weighted aggregation requires nonempty {numerator_col} and {denominator_col} in {summary_path}"
                    )
                numerator = sum(float(row[numerator_col]) for row in group)
                denominator = sum(float(row[denominator_col]) for row in group)
                value = numerator / denominator if denominator else math.nan
            else:
                value = sum(float(row["metric_value"]) for row in group) / len(group)
            units.append({
                "analysis_unit": species, "biological_species": species, "ploidy": ploidy,
                "metric": metric, "metric_value": value, "n_haplotypes_or_subgenomes": len(group), "aggregation": aggregation,
            })
    group_a = [row for row in units if str(row["ploidy"]) in set(polyploid_labels)]
    group_b = [row for row in units if str(row["ploidy"]) == diploid_label]
    if not group_a or not group_b:
        raise SchemaError(
            f"need at least one polyploid ({', '.join(polyploid_labels)}) and one {diploid_label} analysis unit; "
            f"found ploidies {sorted({str(row['ploidy']) for row in units})}"
        )
    values_a = [float(row["metric_value"]) for row in group_a]
    values_b = [float(row["metric_value"]) for row in group_b]
    result = stats.mannwhitneyu(values_a, values_b, alternative="two-sided", method="auto")
    # SciPy reports U for the first group. Rank-biserial sign is explicit below.
    rank_biserial = 1 - 2 * float(result.statistic) / (len(values_a) * len(values_b))
    test_rows = [{
        "test": "Mann-Whitney U (two-sided)", "metric": metric, "analysis_unit": unit,
        "aggregation": aggregation if unit == "species" else "none_haplotype_sensitivity",
        "group_a": "+".join(polyploid_labels), "group_b": diploid_label,
        "n_group_a": len(values_a), "n_group_b": len(values_b), "U_group_a": result.statistic,
        "p_value": result.pvalue, "rank_biserial_group_a_minus_group_b": rank_biserial,
        "interpretation_warning": (
            "Haplotype-level result is a sensitivity analysis; related subgenomes are not independent biological species."
            if unit == "haplotype" else ""
        ),
    }]
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    units_path = out / "ploidy_analysis_units.tsv"
    test_path = out / "ploidy_mann_whitney.tsv"
    write_tsv(units_path, units, ["analysis_unit", "biological_species", "ploidy", "metric", "metric_value", "n_haplotypes_or_subgenomes", "aggregation"])
    write_tsv(test_path, test_rows, ["test", "metric", "analysis_unit", "aggregation", "group_a", "group_b", "n_group_a", "n_group_b", "U_group_a", "p_value", "rank_biserial_group_a_minus_group_b", "interpretation_warning"])
    return {"analysis_units": units_path, "test": test_path}


def _shapiro_normal(values: list[float], stats) -> tuple[bool, float | None]:
    if len(values) < 3 or len(values) > 5000:
        return False, None
    result = stats.shapiro(values)
    return bool(result.pvalue > 0.05), float(result.pvalue)


def _cohen_dz(a: list[float], b: list[float]) -> float | None:
    differences = [x - y for x, y in zip(a, b)]
    if len(differences) < 2:
        return None
    mean = sum(differences) / len(differences)
    variance = sum((value - mean) ** 2 for value in differences) / (len(differences) - 1)
    return mean / math.sqrt(variance) if variance > 0 else None


def _holm_adjust(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    count = len(indexed)
    result = [0.0] * count
    running = 0.0
    for rank, (index, value) in enumerate(indexed, start=1):
        candidate = min(1.0, (count - rank + 1) * value)
        running = max(running, candidate)
        result[index] = running
    return result


def subgenome_comparison(
    spatial_inter_path: str | Path,
    metadata_path: str | Path,
    output_dir: str | Path,
    metric: str = "loss_fragment_per_target_gene",
    method: str = "auto",
) -> dict[str, Path]:
    """Paired/rm comparison across subgenomes using chromosomes as paired units.

    The input is the output from ``spatial-summary``.  The method never
    fabricates a pair: a chromosome must be present in every compared
    subgenome, and excluded chromosomes are written to the matrix file.
    """
    if method not in {"auto", "paired-t", "wilcoxon", "rm-anova", "friedman"}:
        raise SchemaError("method must be auto, paired-t, wilcoxon, rm-anova, or friedman")
    stats = _scipy_stats()
    rows = read_tsv(spatial_inter_path, required=["sample_id", "target_chromosome", metric])
    metadata = _metadata_by_sample(metadata_path)
    label_column = "haplotype_or_subgenome"
    merged: list[dict[str, object]] = []
    for row in rows:
        sample = row["sample_id"]
        meta = metadata.get(sample)
        if meta is None:
            raise SchemaError(f"{spatial_inter_path}: sample_id {sample!r} absent from metadata")
        label = meta.get(label_column) or meta.get("subgenome") or meta.get("haplotype") or ""
        if not label or label.upper() == "NA":
            continue
        if not row[metric]:
            continue
        merged.append({
            "biological_species": meta["biological_species"], "sample_id": sample, "subgenome": label,
            "target_chromosome": row["target_chromosome"], "metric_value": float(row[metric]),
        })
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in merged:
        grouped[str(row["biological_species"])].append(row)
    omnibus: list[dict[str, object]] = []
    pairwise: list[dict[str, object]] = []
    matrices: list[dict[str, object]] = []
    for species, group in sorted(grouped.items()):
        label_to_sample: dict[str, str] = {}
        matrix: dict[str, dict[str, float]] = defaultdict(dict)
        for row in group:
            sub = str(row["subgenome"])
            previous = matrix[str(row["target_chromosome"])].get(sub)
            if previous is not None:
                raise SchemaError(f"{species}: duplicate ({row['target_chromosome']}, {sub}) observations")
            matrix[str(row["target_chromosome"])][sub] = float(row["metric_value"])
            label_to_sample[sub] = str(row["sample_id"])
        labels = sorted(label_to_sample)
        if len(labels) < 2:
            continue
        complete_chromosomes = [chrom for chrom, values in matrix.items() if all(label in values for label in labels)]
        complete_chromosomes.sort(key=lambda value: value)
        for chromosome, values in matrix.items():
            matrices.append({
                "biological_species": species, "target_chromosome": chromosome,
                "is_complete_pairing": str(chromosome in complete_chromosomes).lower(),
                **{f"subgenome_{label}": values.get(label, "") for label in labels},
            })
        n = len(complete_chromosomes)
        if n < 2:
            omnibus.append({"biological_species": species, "test": "not_run", "reason": "fewer_than_two_fully_paired_chromosomes", "n_subgenomes": len(labels), "n_paired_chromosomes": n})
            continue
        vectors = {label: [matrix[chromosome][label] for chromosome in complete_chromosomes] for label in labels}
        if len(labels) == 2:
            a, b = labels
            normal, shapiro_p = _shapiro_normal([x - y for x, y in zip(vectors[a], vectors[b])], stats)
            selected = method
            if selected == "auto":
                selected = "paired-t" if normal else "wilcoxon"
            if selected not in {"paired-t", "wilcoxon"}:
                raise SchemaError(f"{species}: two subgenomes require paired-t or wilcoxon, not {selected}")
            if selected == "paired-t":
                result = stats.ttest_rel(vectors[a], vectors[b])
                statistic, p_value = float(result.statistic), float(result.pvalue)
                effect = _cohen_dz(vectors[a], vectors[b])
                effect_name = "Cohen_dz"
                test_name = "paired_t"
            else:
                try:
                    result = stats.wilcoxon(vectors[a], vectors[b], zero_method="wilcox")
                    statistic, p_value = float(result.statistic), float(result.pvalue)
                except ValueError:  # all differences zero
                    statistic, p_value = 0.0, 1.0
                effect, effect_name, test_name = None, "not_computed", "Wilcoxon_signed_rank"
            omnibus.append({
                "biological_species": species, "test": test_name, "reason": "", "n_subgenomes": 2,
                "n_paired_chromosomes": n, "subgenomes": f"{a} vs {b}", "statistic": statistic,
                "p_value": p_value, "effect_name": effect_name, "effect": effect if effect is not None else "",
                "shapiro_p_of_paired_differences": shapiro_p if shapiro_p is not None else "",
            })
            pairwise.append({
                "biological_species": species, "subgenome_a": a, "subgenome_b": b, "test": test_name,
                "n_paired_chromosomes": n, "statistic": statistic, "p_value_raw": p_value, "p_value_holm": p_value,
                "effect_name": effect_name, "effect": effect if effect is not None else "",
            })
        else:
            all_normal = all(_shapiro_normal([x - y for x, y in zip(vectors[a], vectors[b])], stats)[0] for a, b in itertools.combinations(labels, 2))
            selected = method
            if selected == "auto":
                selected = "rm-anova" if all_normal else "friedman"
            if selected == "rm-anova":
                try:
                    import pandas as pd
                    from statsmodels.stats.anova import AnovaRM
                except ImportError as exc:
                    raise RuntimeError("RM-ANOVA needs statsmodels and pandas; use --method friedman or install statistics extras") from exc
                long = pd.DataFrame([
                    {"chromosome": chromosome, "subgenome": label, "metric": matrix[chromosome][label]}
                    for chromosome in complete_chromosomes for label in labels
                ])
                result_table = AnovaRM(long, depvar="metric", subject="chromosome", within=["subgenome"]).fit().anova_table
                row = result_table.loc["subgenome"]
                statistic, p_value = float(row["F Value"]), float(row["Pr > F"])
                test_name = "RM_ANOVA"
            elif selected == "friedman":
                result = stats.friedmanchisquare(*(vectors[label] for label in labels))
                statistic, p_value, test_name = float(result.statistic), float(result.pvalue), "Friedman"
            else:
                raise SchemaError(f"{species}: >=3 subgenomes require rm-anova or friedman, not {selected}")
            omnibus.append({
                "biological_species": species, "test": test_name, "reason": "", "n_subgenomes": len(labels),
                "n_paired_chromosomes": n, "subgenomes": "+".join(labels), "statistic": statistic, "p_value": p_value,
                "effect_name": "", "effect": "", "shapiro_p_of_paired_differences": "",
            })
            if p_value < 0.05:
                raw_rows: list[dict[str, object]] = []
                raw_p: list[float] = []
                for a, b in itertools.combinations(labels, 2):
                    if test_name == "RM_ANOVA":
                        result = stats.ttest_rel(vectors[a], vectors[b])
                        test = "paired_t"
                        effect, effect_name = _cohen_dz(vectors[a], vectors[b]), "Cohen_dz"
                    else:
                        try:
                            result = stats.wilcoxon(vectors[a], vectors[b], zero_method="wilcox")
                            test = "Wilcoxon_signed_rank"
                        except ValueError:
                            class Result: statistic, pvalue = 0.0, 1.0
                            result, test = Result(), "Wilcoxon_signed_rank"
                        effect, effect_name = None, "not_computed"
                    raw_rows.append({
                        "biological_species": species, "subgenome_a": a, "subgenome_b": b, "test": test,
                        "n_paired_chromosomes": n, "statistic": float(result.statistic), "p_value_raw": float(result.pvalue),
                        "effect_name": effect_name, "effect": effect if effect is not None else "",
                    })
                    raw_p.append(float(result.pvalue))
                corrected = _holm_adjust(raw_p)
                for row, p_corrected in zip(raw_rows, corrected):
                    row["p_value_holm"] = p_corrected
                pairwise.extend(raw_rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    omnibus_path = out / "subgenome_omnibus.tsv"
    pairwise_path = out / "subgenome_pairwise.tsv"
    matrix_path = out / "subgenome_paired_chromosome_matrix.tsv"
    write_tsv(omnibus_path, omnibus, ["biological_species", "test", "reason", "n_subgenomes", "n_paired_chromosomes", "subgenomes", "statistic", "p_value", "effect_name", "effect", "shapiro_p_of_paired_differences"])
    write_tsv(pairwise_path, pairwise, ["biological_species", "subgenome_a", "subgenome_b", "test", "n_paired_chromosomes", "statistic", "p_value_raw", "p_value_holm", "effect_name", "effect"])
    matrix_fields = ["biological_species", "target_chromosome", "is_complete_pairing"] + sorted({key for row in matrices for key in row if key.startswith("subgenome_")})
    write_tsv(matrix_path, matrices, matrix_fields)
    return {"omnibus": omnibus_path, "pairwise": pairwise_path, "paired_matrix": matrix_path}
