#!/usr/bin/env python3
"""Prepare full-unit inputs for all-decayed expression and copy analyses.

Shared and non-shared genes are deliberately pooled.  ``decayed`` is the only
positive class; ``retained``, ``decayed``, and ``deleted`` remain in the
resolved opportunity denominator, while ``not_called_loss`` is excluded.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


CLASSES = {"retained", "decayed", "deleted", "not_called_loss"}
RESOLVED = {"retained", "decayed", "deleted"}


class InputError(ValueError):
    """Raised when the complete matrix and supporting ledgers do not close."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-matrix", required=True, type=Path)
    parser.add_argument("--unit-ledger", required=True, type=Path)
    parser.add_argument("--reference-expression", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise InputError(f"missing or empty input: {path}")
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def read_tsv(path: Path) -> pd.DataFrame:
    binding(path)
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        compression="infer",
    )


def require(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise InputError(f"{label} missing columns: {', '.join(missing)}")


def write_deterministic_gzip(frame: pd.DataFrame, path: Path) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                frame.to_csv(handle, sep="\t", index=False, lineterminator="\n")


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise InputError(f"refusing to overwrite output directory: {args.output_dir}")

    matrix = read_tsv(args.manuscript_matrix)
    require(
        matrix,
        {
            "reference_gene_id",
            "assembly_unit_id",
            "manuscript_classification",
        },
        "manuscript matrix",
    )
    matrix = matrix[
        ["reference_gene_id", "assembly_unit_id", "manuscript_classification"]
    ].copy()
    expected_rows = args.expected_units * args.expected_reference_genes
    if len(matrix) != expected_rows:
        raise InputError(
            f"matrix has {len(matrix)} rows; expected {expected_rows}"
        )
    if matrix.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
        raise InputError("matrix contains duplicate gene-unit rows")
    observed_classes = set(matrix["manuscript_classification"])
    if not observed_classes.issubset(CLASSES):
        raise InputError(f"unsupported classes: {sorted(observed_classes - CLASSES)}")
    reference_genes = set(matrix["reference_gene_id"])
    units = set(matrix["assembly_unit_id"])
    if (
        len(reference_genes) != args.expected_reference_genes
        or len(units) != args.expected_units
    ):
        raise InputError("matrix gene or unit universe does not close")

    ledger = read_tsv(args.unit_ledger)
    require(
        ledger,
        {"sample_id", "ploidy", "analysis_role", "input_scope"},
        "unit ledger",
    )
    ledger = ledger.loc[ledger["analysis_role"] == "target_repertoire"].copy()
    if (
        len(ledger) != args.expected_units
        or ledger["sample_id"].duplicated().any()
        or set(ledger["sample_id"]) != units
    ):
        raise InputError("target unit ledger does not match the matrix")
    if set(ledger["input_scope"]) != {"whole_genome"}:
        raise InputError("all target units must use whole_genome input scope")
    ploidy = dict(zip(ledger["sample_id"], ledger["ploidy"]))

    expression = read_tsv(args.reference_expression)
    require(expression, {"reference_gene_id", "leaf_raw_count"}, "expression")
    expression = expression[["reference_gene_id", "leaf_raw_count"]].copy()
    if (
        expression["reference_gene_id"].duplicated().any()
        or set(expression["reference_gene_id"]) != reference_genes
    ):
        raise InputError("expression gene universe does not match the matrix")
    counts = pd.to_numeric(expression["leaf_raw_count"], errors="coerce")
    if counts.isna().any() or (counts < 0).any():
        raise InputError("leaf_raw_count must be complete and non-negative")

    resolved = matrix.loc[
        matrix["manuscript_classification"].isin(RESOLVED)
    ].copy()
    resolved["ploidy"] = resolved["assembly_unit_id"].map(ploidy)
    resolved = resolved.rename(
        columns={"manuscript_classification": "classification"}
    )
    resolved = resolved[
        ["reference_gene_id", "assembly_unit_id", "ploidy", "classification"]
    ]
    class_counts = Counter(matrix["manuscript_classification"])
    if len(resolved) != expected_rows - class_counts["not_called_loss"]:
        raise InputError("resolved denominator does not close")

    parent = args.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=parent)
    )
    try:
        loss_path = temporary / "resolved_all_unit_loss_table.tsv.gz"
        expression_path = temporary / "all_reference_leaf_raw_counts.tsv"
        write_deterministic_gzip(resolved, loss_path)
        expression.to_csv(
            expression_path,
            sep="\t",
            index=False,
            lineterminator="\n",
        )
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_ALL_DECAYED_EXPRESSION_COPY_INPUTS",
            "created_at_utc": datetime.now(timezone.utc).isoformat(
                timespec="seconds"
            ),
            "definition": {
                "scope": "shared_and_nonshared_pooled",
                "positive": "decayed_only",
                "resolved_denominator": "retained + decayed + deleted",
                "excluded": "not_called_loss",
                "analysis_unit": "independent_assembly_unit",
            },
            "counts": {
                "reference_genes": len(reference_genes),
                "assembly_units": len(units),
                "complete_matrix_rows": len(matrix),
                "resolved_rows": len(resolved),
                "not_called_rows": class_counts["not_called_loss"],
                "retained_rows": class_counts["retained"],
                "decayed_rows": class_counts["decayed"],
                "deleted_rows": class_counts["deleted"],
            },
            "inputs": [
                binding(args.manuscript_matrix),
                binding(args.unit_ledger),
                binding(args.reference_expression),
            ],
            "outputs": [
                binding(loss_path),
                binding(expression_path),
            ],
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (OSError, InputError, pd.errors.ParserError) as error:
        raise SystemExit(f"ERROR: {error}")
    print(f"Wrote all-decayed expression/copy inputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
