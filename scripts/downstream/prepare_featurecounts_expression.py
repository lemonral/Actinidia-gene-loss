#!/usr/bin/env python3
"""Prepare one featureCounts sample on an explicit reference-gene universe.

The ``legacy_featurecounts_fpkm`` mode repairs the archived workflow without
pretending that it is a new transcript-quantification run.  It calculates

    count * 1e9 / (featureCounts Length * total assigned counts)

where ``Length`` is the span reported for ``-t mRNA -g ID``.  That span is not
an exon-summed transcript length, so the result is deliberately labelled as a
legacy intended FPKM sensitivity measurement.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_reference_ids(path: Path) -> list[str]:
    values: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            value = line.rstrip("\n").split("\t", 1)[0].split(None, 1)[0]
            if value.lower() in {"transcript_id", "reference_gene_id", "gene_id", "id"}:
                continue
            values.append(value)
    if not values:
        raise ValueError(f"{path}: no reference identifiers")
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate reference identifiers")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--reference-coords", required=True, type=Path)
    parser.add_argument("--gene-column", default="Geneid")
    parser.add_argument("--length-column", default="Length")
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--count-column", help="Exact featureCounts column header")
    selector.add_argument("--count-column-contains", help="Text that must match exactly one column header")
    parser.add_argument(
        "--measurement",
        choices=["raw_count", "legacy_featurecounts_fpkm"],
        default="raw_count",
        help="FPKM mode requires an exact count header and the matching featureCounts summary",
    )
    parser.add_argument("--summary", type=Path, help="Matching gene_counts.txt.summary (required for FPKM mode)")
    parser.add_argument("--summary-status", default="Assigned")
    parser.add_argument("--output-column", help="Default depends on --measurement")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def require_nonnegative_integer(value: str, label: str) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"{label} must be a non-negative integer, not {value!r}")
    return int(value)


def require_positive_integer(value: str, label: str) -> int:
    parsed = require_nonnegative_integer(value, label)
    if parsed <= 0:
        raise ValueError(f"{label} must be positive, not {value!r}")
    return parsed


def read_summary_assigned(path: Path, count_column: str, status: str) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError(f"{path}: duplicate summary headers")
        if "Status" not in fieldnames:
            raise ValueError(f"{path}: featureCounts summary lacks Status column")
        if count_column not in fieldnames:
            raise ValueError(f"{path}: summary lacks the exact selected sample header")
        matches = [row for row in reader if row.get("Status", "").strip() == status]
    if len(matches) != 1:
        raise ValueError(f"{path}: expected exactly one {status!r} row; found {len(matches)}")
    return require_nonnegative_integer(
        matches[0].get(count_column, "").strip(),
        f"{path}: {status} count for selected sample",
    )


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    if args.measurement == "legacy_featurecounts_fpkm":
        if not args.count_column:
            raise SystemExit("ERROR: FPKM mode requires --count-column with the exact featureCounts header")
        if args.summary is None:
            raise SystemExit("ERROR: FPKM mode requires --summary to verify the assigned-count denominator")
    output_column = args.output_column or (
        "leaf_legacy_intended_fpkm"
        if args.measurement == "legacy_featurecounts_fpkm"
        else "leaf_raw_count"
    )
    try:
        reference_ids = read_reference_ids(args.reference_coords)
        with args.counts.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
            fieldnames = reader.fieldnames or []
            if not fieldnames:
                raise ValueError("featureCounts table has no header")
            if len(fieldnames) != len(set(fieldnames)):
                raise ValueError("featureCounts table contains duplicate headers")
            if args.gene_column not in fieldnames:
                raise ValueError(f"featureCounts table lacks gene column: {args.gene_column}")
            if args.measurement == "legacy_featurecounts_fpkm" and args.length_column not in fieldnames:
                raise ValueError(f"featureCounts table lacks length column: {args.length_column}")
            if args.count_column:
                matches = [name for name in fieldnames if name == args.count_column]
            else:
                matches = [name for name in fieldnames if args.count_column_contains in name]
            if len(matches) != 1:
                raise ValueError(f"count-column selector matched {len(matches)} headers; expected exactly one")
            count_column = matches[0]
            counts: dict[str, int] = {}
            lengths: dict[str, int] = {}
            for row in reader:
                gene = row[args.gene_column].strip()
                if not gene:
                    raise ValueError("featureCounts table contains an empty gene ID")
                if gene in counts:
                    raise ValueError(f"featureCounts table contains duplicate gene ID: {gene}")
                value = row[count_column].strip()
                counts[gene] = require_nonnegative_integer(value, f"count for {gene}")
                if args.measurement == "legacy_featurecounts_fpkm":
                    lengths[gene] = require_positive_integer(
                        row[args.length_column].strip(), f"{args.length_column} for {gene}"
                    )
        if not counts:
            raise ValueError("featureCounts table contains no data rows")
        missing = set(reference_ids).difference(counts)
        if missing:
            examples = ", ".join(sorted(missing)[:5])
            raise ValueError(f"featureCounts table lacks {len(missing)} reference IDs ({examples})")
        total_assigned = sum(counts.values())
        summary_assigned: int | None = None
        if args.measurement == "legacy_featurecounts_fpkm":
            if total_assigned <= 0:
                raise ValueError("selected sample has a zero assigned-count denominator")
            summary_assigned = read_summary_assigned(args.summary, count_column, args.summary_status)
            if summary_assigned != total_assigned:
                raise ValueError(
                    "selected count-column sum does not equal the featureCounts summary denominator: "
                    f"{total_assigned} versus {summary_assigned}"
                )
    except (OSError, ValueError, csv.Error) as exc:
        raise SystemExit(f"ERROR: {exc}")

    reference_set = set(reference_ids)
    extra_ids = set(counts).difference(reference_set)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    if args.measurement == "legacy_featurecounts_fpkm":
        output = args.output_dir / "csc_leaf_legacy_intended_fpkm_canonical.tsv"
        prepared_values = {
            gene: format(counts[gene] * 1_000_000_000 / (lengths[gene] * total_assigned), ".12g")
            for gene in reference_ids
        }
    else:
        output = args.output_dir / "csc_leaf_raw_counts_canonical.tsv"
        prepared_values = {gene: str(counts[gene]) for gene in reference_ids}
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reference_gene_id", output_column])
        writer.writerows((gene, prepared_values[gene]) for gene in reference_ids)
    with (args.output_dir / "preparation_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "value"])
        qc_rows: list[tuple[str, object]] = [
            ("featurecounts_gene_ids", len(counts)),
            ("reference_gene_ids", len(reference_ids)),
            ("reference_ids_missing_count", 0),
            ("noncanonical_featurecounts_ids_excluded", len(extra_ids)),
            ("selected_count_column", count_column),
            ("measurement", args.measurement),
            ("assigned_count_denominator_all_featurecounts_ids", total_assigned),
            ("zero_count_featurecounts_ids", sum(value == 0 for value in counts.values())),
        ]
        if args.measurement == "legacy_featurecounts_fpkm":
            qc_rows.extend([
                ("length_column", args.length_column),
                ("minimum_mrna_genomic_span", min(lengths.values())),
                ("maximum_mrna_genomic_span", max(lengths.values())),
                ("summary_status", args.summary_status),
                ("summary_assigned_count", summary_assigned),
                ("count_sum_equals_summary_assigned", str(total_assigned == summary_assigned).upper()),
            ])
        writer.writerows(qc_rows)
    if args.measurement == "legacy_featurecounts_fpkm":
        measurement_description = (
            "legacy intended FPKM sensitivity: raw_count*1e9/(featureCounts Length*total assigned counts); "
            "Length is the genomic span from featureCounts -t mRNA -g ID, not exon-summed transcript length; "
            "not a Cufflinks-requantified FPKM"
        )
    else:
        measurement_description = "featureCounts raw count; not FPKM"
    metadata = [
        ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
        ("command", " ".join(shlex.quote(value) for value in sys.argv)),
        ("dataset_id", args.dataset_id),
        ("measurement", measurement_description),
        ("selected_count_column", count_column),
        ("assigned_count_denominator", total_assigned),
        ("denominator_scope", "sum across every featureCounts feature ID before canonical-reference filtering"),
        ("counts_sha256", sha256(args.counts)),
        ("summary_sha256", sha256(args.summary) if args.summary else ""),
        ("reference_coords_sha256", sha256(args.reference_coords)),
        ("output_sha256", sha256(output)),
    ]
    with (args.output_dir / "run_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerows(metadata)
    print(f"Wrote {len(reference_ids)} canonical {args.measurement} rows to {output}")


if __name__ == "__main__":
    main()
