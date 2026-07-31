#!/usr/bin/env python3
"""Calculate four-tissue TPM and the arithmetic mean on a frozen gene universe."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


TISSUE_TOKENS = {
    "stem": "S23I0030",
    "leaf": "S23I0033",
    "fruit": "S23I0032",
    "root": "S23I0031",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_attributes(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    return values


def read_reference_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "reference_gene_id" not in (reader.fieldnames or []):
            raise ValueError("reference table lacks reference_gene_id")
        values = [row["reference_gene_id"].strip() for row in reader]
    if not values or any(not value for value in values):
        raise ValueError("reference table contains no IDs or an empty ID")
    if len(values) != len(set(values)):
        raise ValueError("reference table contains duplicate IDs")
    return values


def read_exon_lengths(path: Path, reference_ids: set[str]) -> dict[str, int]:
    intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = parse_attributes(fields[8])
            for parent in attrs.get("Parent", "").split(","):
                if parent in reference_ids:
                    intervals[parent].append((fields[0], int(fields[3]), int(fields[4])))
    lengths: dict[str, int] = {}
    for gene in reference_ids:
        by_chromosome: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chromosome, start, end in intervals.get(gene, []):
            by_chromosome[chromosome].append((start, end))
        total = 0
        for chromosome_intervals in by_chromosome.values():
            merged: list[list[int]] = []
            for start, end in sorted(chromosome_intervals):
                if not merged or start > merged[-1][1] + 1:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            total += sum(end - start + 1 for start, end in merged)
        if total <= 0:
            raise ValueError(f"no positive exon-union length for {gene}")
        lengths[gene] = total
    return lengths


def read_counts(path: Path, reference_ids: set[str]) -> tuple[dict[str, dict[str, int]], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        headers = reader.fieldnames or []
        selected: dict[str, str] = {}
        for tissue, token in TISSUE_TOKENS.items():
            matches = [header for header in headers if token in header]
            if len(matches) != 1:
                raise ValueError(f"{tissue}: expected one count column containing {token}; found {len(matches)}")
            selected[tissue] = matches[0]
        values: dict[str, dict[str, int]] = {}
        for row in reader:
            gene = row.get("Geneid", "").strip()
            if gene not in reference_ids:
                continue
            if gene in values:
                raise ValueError(f"duplicate count row for {gene}")
            values[gene] = {tissue: int(row[column]) for tissue, column in selected.items()}
    missing = reference_ids.difference(values)
    if missing:
        raise ValueError(f"counts lack {len(missing)} reference IDs; first={sorted(missing)[0]}")
    return values, selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--reference-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: output directory already exists: {args.output_dir}")
    try:
        reference_ids = read_reference_ids(args.reference_table)
        reference_set = set(reference_ids)
        exon_lengths = read_exon_lengths(args.gff, reference_set)
        counts, selected_columns = read_counts(args.counts, reference_set)
        rpk = {
            tissue: {
                gene: counts[gene][tissue] / (exon_lengths[gene] / 1000.0)
                for gene in reference_ids
            }
            for tissue in TISSUE_TOKENS
        }
        scales = {tissue: sum(values.values()) / 1_000_000.0 for tissue, values in rpk.items()}
        if any(scale <= 0 for scale in scales.values()):
            raise ValueError("one or more tissue RPK denominators are zero")
        tpm = {
            tissue: {gene: rpk[tissue][gene] / scales[tissue] for gene in reference_ids}
            for tissue in TISSUE_TOKENS
        }
    except (OSError, ValueError, csv.Error) as exc:
        raise SystemExit(f"ERROR: {exc}")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = args.output_dir / "csc_four_tissue_mean_tpm.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "reference_gene_id",
            "stem_tpm",
            "leaf_tpm",
            "fruit_tpm",
            "root_tpm",
            "four_tissue_mean_tpm",
        ])
        for gene in reference_ids:
            values = [tpm[tissue][gene] for tissue in ("stem", "leaf", "fruit", "root")]
            writer.writerow([gene, *(format(value, ".12g") for value in values), format(sum(values) / 4.0, ".12g")])

    qc = args.output_dir / "preparation_qc.tsv"
    with qc.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "value"])
        writer.writerow(["reference_gene_count", len(reference_ids)])
        writer.writerow(["minimum_exon_union_length_bp", min(exon_lengths.values())])
        writer.writerow(["maximum_exon_union_length_bp", max(exon_lengths.values())])
        for tissue in ("stem", "leaf", "fruit", "root"):
            writer.writerow([f"{tissue}_count_column", selected_columns[tissue]])
            writer.writerow([f"{tissue}_tpm_sum", format(sum(tpm[tissue].values()), ".12g")])

    metadata = args.output_dir / "run_metadata.tsv"
    with metadata.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerows([
            ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
            ("measurement", "TPM from archived featureCounts counts and transcript exon-union lengths"),
            ("cross_tissue_summary", "arithmetic mean of stem, leaf, fruit, and root TPM"),
            ("biological_replicates_per_tissue", "1"),
            ("counts_sha256", sha256(args.counts)),
            ("gff_sha256", sha256(args.gff)),
            ("reference_table_sha256", sha256(args.reference_table)),
            ("output_sha256", sha256(output)),
        ])
    print(f"Wrote {len(reference_ids)} four-tissue TPM rows to {output}")


if __name__ == "__main__":
    main()
