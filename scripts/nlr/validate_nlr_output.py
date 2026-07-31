#!/usr/bin/env python3
"""Validate NLR-Annotator outputs and produce raw/major-class count tables.

Use the raw eight NLR-Annotator architecture labels as the primary audit
output. A separate class crosswalk may produce a manuscript-facing major-class
summary, but it cannot manufacture RNL or NLR-ID evidence that is absent from
the raw/domain annotation inputs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


DEFAULT_RAW_CLASSES = ",".join([
    "CC-NBARC", "CC-NBARC-LRR", "NBARC", "NBARC-LRR",
    "TIR", "TIR-LRR", "TIR-NBARC", "TIR-NBARC-LRR",
])


def read_crosswalk(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = ["raw_architecture", "major_class", "complete_nbarc_lrr"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Crosswalk missing: {', '.join(missing)}")
    if frame["raw_architecture"].duplicated().any():
        raise ValueError("Crosswalk has duplicate raw_architecture rows")
    return frame[required]


def sample_name(path: Path, suffix: str) -> str:
    name = path.name
    return name[:-len(suffix)] if suffix and name.endswith(suffix) else path.stem


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True,
                        help="NLR-Annotator tab-delimited output; repeat once per haplotype")
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--sample-suffix", default="_output.txt")
    parser.add_argument("--expected-samples", type=int, default=21)
    parser.add_argument("--expected-raw-classes", default=DEFAULT_RAW_CLASSES)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    try:
        crosswalk = read_crosswalk(args.crosswalk)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"ERROR: {exc}")
    mapping = crosswalk.set_index("raw_architecture")["major_class"].to_dict()
    expected = {item.strip() for item in args.expected_raw_classes.split(",") if item.strip()}
    parsed: list[pd.DataFrame] = []
    qc: list[dict[str, object]] = []
    errors: list[str] = []

    for path in args.input:
        try:
            frame = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
        except (OSError, pd.errors.ParserError) as exc:
            errors.append(f"Could not read {path}: {exc}")
            continue
        label = sample_name(path, args.sample_suffix)
        if frame.shape[1] < 3:
            errors.append(f"{path} has fewer than three columns")
            continue
        frame = frame.iloc[:, :3].copy()
        frame.columns = ["gene_id", "nlr_locus_id", "raw_architecture"]
        for column in frame.columns:
            frame[column] = frame[column].astype(str).str.strip()
        duplicates = int(frame["gene_id"].duplicated(keep=False).sum())
        empty = int((frame["gene_id"] == "").sum() + (frame["raw_architecture"] == "").sum())
        unknown = sorted(set(frame["raw_architecture"]) - set(mapping))
        if duplicates:
            errors.append(f"{label} has {duplicates} duplicate GeneID values")
        if empty:
            errors.append(f"{label} has {empty} empty gene/class fields")
        if unknown:
            errors.append(f"{label} has unmapped architectures: {', '.join(unknown)}")
        frame["sample_id"] = label
        frame["major_class"] = frame["raw_architecture"].map(mapping).fillna("UNMAPPED")
        parsed.append(frame)
        qc.append({
            "sample_id": label,
            "records": len(frame),
            "unique_gene_ids": frame["gene_id"].nunique(),
            "duplicate_gene_id_rows": duplicates,
            "empty_gene_or_class_fields": empty,
            "unmapped_raw_classes": ";".join(unknown),
        })

    labels = [sample_name(path, args.sample_suffix) for path in args.input]
    if len(set(labels)) != len(labels):
        errors.append("Input file names resolve to duplicate sample IDs; adjust --sample-suffix or rename inputs")
    if args.expected_samples is not None and len(labels) != args.expected_samples:
        errors.append(f"Observed {len(labels)} inputs; expected {args.expected_samples}")
    if not parsed:
        errors.append("No valid NLR-Annotator output was parsed")

    combined = pd.concat(parsed, ignore_index=True) if parsed else pd.DataFrame(
        columns=["gene_id", "nlr_locus_id", "raw_architecture", "sample_id", "major_class"]
    )
    observed = set(combined["raw_architecture"])
    extras = sorted(observed - expected)
    if extras:
        errors.append(f"Observed raw architectures outside the declared NLR-Annotator contract: {', '.join(extras)}")
    qc.append({
        "sample_id": "ALL",
        "records": len(combined),
        "unique_gene_ids": combined[["sample_id", "gene_id"]].drop_duplicates().shape[0],
        "duplicate_gene_id_rows": "",
        "empty_gene_or_class_fields": "",
        "unmapped_raw_classes": ";".join(extras),
    })

    args.output_dir.mkdir(parents=True, exist_ok=False)
    combined.to_csv(args.output_dir / "nlr_calls_tidy.tsv", sep="\t", index=False)
    raw_counts = (
        combined.groupby(["sample_id", "raw_architecture"], as_index=False)
        .size().rename(columns={"size": "n_lost_genes"})
    )
    major_counts = (
        combined.groupby(["sample_id", "major_class"], as_index=False)
        .size().rename(columns={"size": "n_lost_genes"})
    )
    raw_counts.to_csv(args.output_dir / "nlr_raw_architecture_counts.tsv", sep="\t", index=False)
    major_counts.to_csv(args.output_dir / "nlr_major_class_counts.tsv", sep="\t", index=False)
    pd.DataFrame(qc).to_csv(args.output_dir / "nlr_qc.tsv", sep="\t", index=False)
    if errors:
        print("FAILED: " + "; ".join(errors), file=sys.stderr)
        raise SystemExit(2)
    print(f"Validated {len(labels)} NLR outputs; results written to {args.output_dir}")


if __name__ == "__main__":
    main()
