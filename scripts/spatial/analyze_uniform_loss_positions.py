#!/usr/bin/env python3
"""Test target-chromosome position of strictly supported pseudogenized genes.

The primary analysis uses only observed target-genome loci: strictly supported
pseudogenized calls versus exact-SynOrths retained calls.  Deleted calls have no
observed target gene and therefore remain an expected-locus sensitivity rather
than entering the primary position numerator or denominator.  Uncertain rows
never enter either.  Unit fixed effects control assembly-unit baselines;
legacy/new interaction and source-stratified fits expose residual
source-dependent behavior despite the identical classification pipeline.
Reference-gene clustered sandwich errors account for repeated observations of
each reference gene across units.
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
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2, norm


class SpatialError(RuntimeError):
    pass


UNIFORM_CONFIG_COLUMNS = ("unit", "target_genome", "candidate_dir", "output_dir")
REQUIRED_MATRIX_COLUMNS = {
    "reference_gene_id", "assembly_unit_id", "source_group", "classification", "callable",
    "positive_loss", "target_chromosome", "position_midpoint_1based",
    "position_start_1based", "position_end_1based", "query_coverage",
    "exact_alignment_identity", "alignment_score", "frameshift_events",
    "inframe_stop_codons", "disruption_supported",
    "resolved_for_spatial_model", "callable_interval_chromosome",
    "callable_interval_midpoint_1based",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.resolve().open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (source.stat().st_size == 0 and not allow_empty):
        raise SpatialError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise SpatialError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise SpatialError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size == 0):
        raise SpatialError(f"missing input file: {value!r}")
    return path


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SpatialError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise SpatialError(f"{path}: JSON root is not an object")
    return value


def checksum_row(directory: Path, filename: str) -> None:
    with (directory / "checksums.tsv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row.get("file") == filename]
    if len(selected) != 1:
        raise SpatialError(f"{directory.name}: checksum row missing for {filename}")
    observed = binding(directory / filename, allow_empty=True)
    if selected[0].get("bytes") != str(observed["bytes"]) or selected[0].get("sha256") != observed["sha256"]:
        raise SpatialError(f"{directory.name}: checksum mismatch for {filename}")


def read_uniform_config(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if tuple(reader.fieldnames or ()) != UNIFORM_CONFIG_COLUMNS or not rows:
        raise SpatialError("uniform config columns/rows differ from exact schema")
    units = [row["unit"] for row in rows]
    if len(units) != len(set(units)) or any(not value for row in rows for value in row.values()):
        raise SpatialError("uniform config contains empty/duplicate values")
    return rows


def read_paf_lengths(path: Path, unit: str) -> dict[str, int]:
    result: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SpatialError(f"{unit}:PAF:{line_number}: fewer than 12 fields")
            chromosome = fields[5]
            try:
                length = int(fields[6])
            except ValueError as error:
                raise SpatialError(f"{unit}:PAF:{line_number}: invalid target length") from error
            if length < 2 or (chromosome in result and result[chromosome] != length):
                raise SpatialError(f"{unit}:PAF:{line_number}: inconsistent target length")
            result[chromosome] = length
    if not result:
        raise SpatialError(f"{unit}: raw PAF has no target-length evidence")
    return result


def normalized_end_distance(midpoint: float, length: int, context: str) -> float:
    if not 1 <= midpoint <= length:
        raise SpatialError(f"{context}: midpoint lies outside chromosome")
    return min(midpoint - 1.0, length - midpoint) / ((length - 1.0) / 2.0)


def design(unit: np.ndarray, x: np.ndarray, source: np.ndarray, *, interaction: bool) -> tuple[np.ndarray, list[str]]:
    units = sorted(set(int(value) for value in unit))
    if units != list(range(len(units))):
        raise SpatialError("unit indices are not contiguous")
    columns = [np.ones(len(x), dtype=float)]
    names = ["intercept"]
    for index in units[1:]:
        columns.append((unit == index).astype(float))
        names.append(f"unit_fixed_effect_{index}")
    columns.append(x)
    names.append("normalized_end_distance")
    if interaction:
        columns.append(x * source)
        names.append("new_source_x_normalized_end_distance")
    return np.column_stack(columns), names


def fit_logistic(
    unit: np.ndarray,
    x: np.ndarray,
    source: np.ndarray,
    y: np.ndarray,
    gene: np.ndarray,
    *,
    interaction: bool,
) -> dict[str, object]:
    matrix, names = design(unit, x, source, interaction=interaction)
    if not 0 < y.sum() < len(y):
        raise SpatialError("logistic response has only one class")
    start = np.zeros(matrix.shape[1], dtype=float)
    start[0] = math.log(float(y.mean()) / (1.0 - float(y.mean())))

    def objective(beta):
        eta = matrix @ beta
        return float(np.logaddexp(0.0, eta).sum() - y @ eta)

    def gradient(beta):
        return matrix.T @ (expit(matrix @ beta) - y)

    fitted = minimize(objective, start, jac=gradient, method="L-BFGS-B", options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-8})
    if not fitted.success or not np.all(np.isfinite(fitted.x)):
        raise SpatialError(f"logistic fit failed: {fitted.message}")
    beta = fitted.x
    probability = expit(matrix @ beta)
    residual = y - probability
    weights = probability * (1.0 - probability)
    information = matrix.T @ (weights[:, None] * matrix)
    bread = np.linalg.pinv(information, rcond=1e-12)
    gene_count = int(gene.max()) + 1
    clustered_scores = np.column_stack(
        [np.bincount(gene, weights=residual * matrix[:, index], minlength=gene_count) for index in range(matrix.shape[1])]
    )
    correction = (gene_count / (gene_count - 1.0)) * ((len(y) - 1.0) / (len(y) - matrix.shape[1]))
    covariance = bread @ (clustered_scores.T @ clustered_scores) @ bread * correction
    standard_error = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    coefficients = []
    for index, name in enumerate(names):
        estimate = beta[index]
        error = standard_error[index]
        z = estimate / error if error > 0 else math.nan
        coefficients.append(
            {
                "term": name,
                "estimate_log_odds": estimate,
                "cluster_robust_se": error,
                "odds_ratio": float(np.exp(np.clip(estimate, -700, 700))),
                "ci95_low_log_odds": estimate - 1.959963984540054 * error,
                "ci95_high_log_odds": estimate + 1.959963984540054 * error,
                "wald_p": 2.0 * norm.sf(abs(z)) if math.isfinite(z) else math.nan,
            }
        )
    return {
        "coefficients": coefficients,
        "negative_log_likelihood": float(fitted.fun),
        "n": len(y),
        "positive": int(y.sum()),
        "genes": gene_count,
        "parameters": matrix.shape[1],
    }


def reindex(values: np.ndarray) -> np.ndarray:
    mapping = {value: index for index, value in enumerate(sorted(set(int(item) for item in values)))}
    return np.asarray([mapping[int(item)] for item in values], dtype=np.int16)


def model_specifications(
    class_array: np.ndarray,
    source_array: np.ndarray,
) -> list[tuple[str, np.ndarray, bool]]:
    """Return the frozen primary and sensitivity position-model cohorts."""
    pseudogenized_or_retained = class_array != "deleted"
    deleted_or_retained = class_array != "pseudogenized"
    return [
        (
            "primary_pseudogenized_unit_fe_source_interaction",
            pseudogenized_or_retained,
            True,
        ),
        ("pooled_pseudogenized_vs_retained_unit_fe", pseudogenized_or_retained, False),
        (
            "legacy_pseudogenized_vs_retained_unit_fe",
            pseudogenized_or_retained & (source_array == 0),
            False,
        ),
        (
            "new_pseudogenized_vs_retained_unit_fe",
            pseudogenized_or_retained & (source_array == 1),
            False,
        ),
        ("deleted_expected_locus_vs_retained_sensitivity", deleted_or_retained, False),
        (
            "combined_positive_expected_locus_sensitivity",
            np.ones(len(class_array), dtype=bool),
            False,
        ),
    ]


def write_tsv(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...] | list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", required=True, type=Path)
    parser.add_argument("--uniform-config", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=10)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        if not 4 <= args.bins <= 50:
            raise SpatialError("bins must be between 4 and 50")
        root = args.data_root.resolve()
        matrix_dir = args.matrix_dir.resolve()
        matrix_path = matrix_dir / "complete_unit_loss_matrix.tsv"
        manifest = strict_json(matrix_dir / "run_manifest.json")
        if manifest.get("status") != "PASS" or manifest.get("workflow") != "uniform_old_new_complete_loss_matrix":
            raise SpatialError("matrix manifest is not the exact uniform PASS")
        checksum_row(matrix_dir, "complete_unit_loss_matrix.tsv")
        config_rows = read_uniform_config(args.uniform_config)
        unit_order = [row["unit"] for row in config_rows]
        unit_index = {unit: index for index, unit in enumerate(unit_order)}

        lengths_by_unit: dict[str, dict[str, int]] = {}
        uniform_audit = []
        for row in config_rows:
            unit = row["unit"]
            directory = resolve(root, row["output_dir"])
            run_manifest = strict_json(directory / "run_manifest.json")
            if run_manifest.get("status") != "PASS" or run_manifest.get("unit") != unit:
                raise SpatialError(f"{unit}: uniform output is not exact PASS")
            checksum_row(directory, "raw_alignments.paf.gz")
            lengths_by_unit[unit] = read_paf_lengths(directory / "raw_alignments.paf.gz", unit)
            uniform_audit.append({"unit": unit, "run_manifest": binding(directory / "run_manifest.json"), "raw_paf": binding(directory / "raw_alignments.paf.gz")})

        units: list[int] = []
        sources: list[int] = []
        genes: list[int] = []
        positions: list[float] = []
        responses: list[int] = []
        classes: list[str] = []
        pseudogenized_positions: list[dict[str, object]] = []
        bins: Counter[tuple[str, str, int, str]] = Counter()
        interval_rows: list[tuple[int, int, int, float, int]] = []
        gene_index: dict[str, int] = {}
        seen_pairs: set[tuple[str, str]] = set()
        with matrix_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not REQUIRED_MATRIX_COLUMNS.issubset(reader.fieldnames or []):
                raise SpatialError("matrix lacks required spatial columns")
            for line_number, row in enumerate(reader, 2):
                unit = row["assembly_unit_id"]
                gene_name = row["reference_gene_id"]
                if unit not in unit_index or row["source_group"] not in {"legacy", "new"}:
                    raise SpatialError(f"matrix:{line_number}: unknown unit/source")
                pair = (unit, gene_name)
                if pair in seen_pairs:
                    raise SpatialError(f"matrix:{line_number}: duplicate unit/gene")
                seen_pairs.add(pair)
                gene_id = gene_index.setdefault(gene_name, len(gene_index))
                source_id = int(row["source_group"] == "new")
                classification = row["classification"]
                if row["resolved_for_spatial_model"] == "true":
                    if classification not in {"retained", "deleted", "pseudogenized"}:
                        raise SpatialError(f"matrix:{line_number}: invalid resolved class")
                    chromosome = row["target_chromosome"]
                    if chromosome not in lengths_by_unit[unit]:
                        raise SpatialError(f"matrix:{line_number}: chromosome length absent for {unit}/{chromosome}")
                    midpoint = float(row["position_midpoint_1based"])
                    x = normalized_end_distance(midpoint, lengths_by_unit[unit][chromosome], f"{unit}/{gene_name}")
                    y = int(classification in {"deleted", "pseudogenized"})
                    units.append(unit_index[unit]); sources.append(source_id); genes.append(gene_id)
                    positions.append(x); responses.append(y); classes.append(classification)
                    bin_number = min(int(x * args.bins), args.bins - 1) + 1
                    bins[(unit, row["source_group"], bin_number, classification)] += 1
                    if classification == "pseudogenized":
                        if row["disruption_supported"] != "true":
                            raise SpatialError(
                                f"matrix:{line_number}: pseudogenized row lacks strict disruption support"
                            )
                        pseudogenized_positions.append(
                            {
                                "reference_gene_id": gene_name,
                                "assembly_unit_id": unit,
                                "source_group": row["source_group"],
                                "target_chromosome": chromosome,
                                "position_start_1based": row["position_start_1based"],
                                "position_end_1based": row["position_end_1based"],
                                "position_midpoint_1based": row["position_midpoint_1based"],
                                "normalized_end_distance": x,
                                "end_to_center_bin": bin_number,
                                "query_coverage": row["query_coverage"],
                                "exact_alignment_identity": row["exact_alignment_identity"],
                                "alignment_score": row["alignment_score"],
                                "frameshift_events": row["frameshift_events"],
                                "inframe_stop_codons": row["inframe_stop_codons"],
                            }
                        )
                if (
                    classification in {"deleted", "pseudogenized", "uncertain"}
                    and row["callable"] == "true"
                    and row["callable_interval_chromosome"]
                    and row["callable_interval_midpoint_1based"]
                ):
                    chromosome = row["callable_interval_chromosome"]
                    if chromosome not in lengths_by_unit[unit]:
                        raise SpatialError(f"matrix:{line_number}: interval chromosome length absent")
                    x = normalized_end_distance(float(row["callable_interval_midpoint_1based"]), lengths_by_unit[unit][chromosome], f"{unit}/{gene_name}:interval")
                    interval_rows.append((unit_index[unit], source_id, gene_id, x, int(classification != "uncertain")))

        expected_pairs = int(manifest["parameters"]["matrix_row_count"])
        if len(seen_pairs) != expected_pairs or set(unit_order) != {pair[0] for pair in seen_pairs}:
            raise SpatialError("matrix pair-count/unit closure failed")
        unit_array = np.asarray(units, dtype=np.int16)
        source_array = np.asarray(sources, dtype=np.int8)
        gene_array = np.asarray(genes, dtype=np.int32)
        x_array = np.asarray(positions, dtype=float)
        y_array = np.asarray(responses, dtype=float)
        class_array = np.asarray(classes)

        model_rows: list[dict[str, object]] = []
        coefficient_rows: list[dict[str, object]] = []
        specifications = model_specifications(class_array, source_array)
        for name, selected, interaction in specifications:
            fit = fit_logistic(
                reindex(unit_array[selected]), x_array[selected], source_array[selected], y_array[selected],
                np.unique(gene_array[selected], return_inverse=True)[1].astype(np.int32), interaction=interaction,
            )
            model_rows.append({
                "model": name, "resolved_rows": fit["n"], "positive_rows": fit["positive"],
                "retained_rows": int(fit["n"]) - int(fit["positive"]), "reference_gene_clusters": fit["genes"],
                "parameters": fit["parameters"], "negative_log_likelihood": fit["negative_log_likelihood"],
                "uncertain_rows_in_model": 0,
            })
            for coefficient in fit["coefficients"]:
                coefficient_rows.append({"model": name, **coefficient})

        interval_array = np.asarray(interval_rows, dtype=float)
        diagnostic = fit_logistic(
            reindex(interval_array[:, 0].astype(np.int16)), interval_array[:, 3], interval_array[:, 1],
            interval_array[:, 4],
            np.unique(interval_array[:, 2].astype(np.int32), return_inverse=True)[1].astype(np.int32),
            interaction=True,
        )
        diagnostic_rows = [{"model": "callable_candidate_positive_vs_uncertain", **row} for row in diagnostic["coefficients"]]

        output = args.output_dir.resolve()
        if output.exists():
            raise SpatialError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        bin_rows = []
        for unit in unit_order:
            source_group = next(row["source_group"] for row in manifest["metrics"]["unit_summaries"] if row["assembly_unit_id"] == unit)
            for bin_number in range(1, args.bins + 1):
                values = {classification: bins[(unit, source_group, bin_number, classification)] for classification in ("retained", "deleted", "pseudogenized")}
                primary_positive = values["pseudogenized"]
                primary_denominator = values["retained"] + primary_positive
                combined_positive = values["deleted"] + primary_positive
                combined_denominator = sum(values.values())
                bin_rows.append({
                    "assembly_unit_id": unit, "source_group": source_group, "end_to_center_bin": bin_number,
                    "normalized_end_distance_lower": (bin_number - 1) / args.bins,
                    "normalized_end_distance_upper": bin_number / args.bins,
                    **values,
                    "primary_pseudogenized": primary_positive,
                    "primary_observed_locus_denominator": primary_denominator,
                    "primary_pseudogenized_rate": (
                        primary_positive / primary_denominator if primary_denominator else ""
                    ),
                    "combined_positive_sensitivity": combined_positive,
                    "combined_expected_locus_denominator": combined_denominator,
                    "combined_positive_sensitivity_rate": (
                        combined_positive / combined_denominator if combined_denominator else ""
                    ),
                })
        length_rows = [
            {"assembly_unit_id": unit, "target_chromosome": chromosome, "chromosome_length_bp": length}
            for unit in unit_order for chromosome, length in sorted(lengths_by_unit[unit].items())
        ]
        write_tsv(staging / "position_bin_summary.tsv", bin_rows, tuple(bin_rows[0]))
        write_tsv(
            staging / "pseudogenized_target_positions.tsv",
            pseudogenized_positions,
            (
                "reference_gene_id", "assembly_unit_id", "source_group", "target_chromosome",
                "position_start_1based", "position_end_1based", "position_midpoint_1based",
                "normalized_end_distance", "end_to_center_bin", "query_coverage",
                "exact_alignment_identity", "alignment_score", "frameshift_events",
                "inframe_stop_codons",
            ),
        )
        write_tsv(staging / "model_summary.tsv", model_rows, tuple(model_rows[0]))
        write_tsv(staging / "model_coefficients.tsv", coefficient_rows, tuple(coefficient_rows[0]))
        write_tsv(staging / "callable_candidate_uncertain_diagnostic.tsv", diagnostic_rows, tuple(diagnostic_rows[0]))
        write_tsv(staging / "chromosome_lengths_from_miniprot_paf.tsv", length_rows, tuple(length_rows[0]))
        run_manifest = {
            "schema_version": 1,
            "workflow": "uniform_pseudogenized_target_position_analysis",
            "status": "PASS",
            "finished_at_utc": utc_now(),
            "definitions": {
                "normalized_end_distance": "0 at chromosome ends; 1 at chromosome center",
                "primary_numerator": "strictly supported pseudogenized calls at observed target-genome alignment loci",
                "primary_denominator": "uniquely positioned retained plus strictly supported pseudogenized observed loci",
                "deleted": "excluded from the primary observed-locus model; expected-locus midpoint used only in labelled sensitivity models",
                "uncertain": "excluded from every positive numerator and from the primary observed-locus denominator",
                "centromere": "not analyzed because no independent centromere intervals are available",
            },
            "statistics": {
                "model": "binomial logistic regression with assembly-unit fixed effects",
                "uncertainty": "reference-gene clustered sandwich standard errors",
                "source_audit": "legacy/new position-slope interaction plus source-stratified fits",
            },
            "inputs": {"matrix_manifest": binding(matrix_dir / "run_manifest.json"), "matrix": binding(matrix_path), "uniform_config": binding(args.uniform_config), "uniform_outputs": uniform_audit},
            "metrics": {
                "matrix_rows": len(seen_pairs),
                "primary_observed_locus_rows": int(np.sum(class_array != "deleted")),
                "primary_pseudogenized_rows": int(np.sum(class_array == "pseudogenized")),
                "primary_retained_rows": int(np.sum(class_array == "retained")),
                "deleted_expected_locus_sensitivity_rows": int(np.sum(class_array == "deleted")),
                "combined_positionable_rows": len(y_array),
                "callable_candidate_diagnostic_rows": len(interval_rows),
            },
        }
        (staging / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (staging / "checksums.tsv").open("w", encoding="utf-8") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output); staging = None
        print(json.dumps({
            "status": "PASS",
            "primary_observed_locus_rows": int(np.sum(class_array != "deleted")),
            "primary_pseudogenized_rows": int(np.sum(class_array == "pseudogenized")),
            "output": str(output),
        }, sort_keys=True))
        return 0
    except (OSError, csv.Error, ValueError, KeyError, SpatialError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    finally:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
