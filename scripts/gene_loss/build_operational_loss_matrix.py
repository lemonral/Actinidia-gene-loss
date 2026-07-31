#!/usr/bin/env python3
"""Materialize a complete operational gene-loss call matrix.

The frozen operational inputs contain positive calls only: one
pseudogenized/decayed list and one deleted list per target haplotype. This
builder expands those lists against an explicit reference-gene universe and
sample manifest. A gene absent from both positive lists is labelled
``not_called_loss``. That label means only that the pipeline made no positive
loss call; it must never be interpreted as proof that the gene is
biologically retained.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from datetime import datetime, timezone
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "t", "yes", "y"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_id_lines(path: Path, *, allow_empty: bool) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not allow_empty and not values:
        raise ValueError(f"{path}: no identifiers")
    duplicates = len(values) - len(set(values))
    if duplicates:
        raise ValueError(f"{path}: {duplicates} duplicate identifier rows")
    if any(len(value.split()) != 1 for value in values):
        raise ValueError(f"{path}: identifiers must be single whitespace-free fields")
    return values


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
    duplicates = len(values) - len(set(values))
    if duplicates:
        raise ValueError(f"{path}: {duplicates} duplicate reference identifiers")
    return values


def read_samples(path: Path, include_column: str, expected_samples: int | None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = {"target_haplotype", "species", "ploidy", include_column}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    selected = [row for row in rows if truthy(row[include_column])]
    if not selected:
        raise ValueError(f"{path}: no rows selected by {include_column}")
    if expected_samples is not None and len(selected) != expected_samples:
        raise ValueError(f"{path}: selected {len(selected)} samples; expected {expected_samples}")
    sample_ids = [row["target_haplotype"].strip() for row in selected]
    if any(not value for value in sample_ids):
        raise ValueError(f"{path}: selected row has an empty target_haplotype")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"{path}: duplicate selected target_haplotype values")
    for row in selected:
        for column in ("target_haplotype", "species", "ploidy"):
            row[column] = row[column].strip()
            if not row[column]:
                raise ValueError(f"{path}: selected sample has an empty {column}")
    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-coords", required=True, type=Path,
                        help="Reference universe; the first field of every data row is the gene ID")
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--include-column", default="include_downstream")
    parser.add_argument("--expected-samples", type=int)
    parser.add_argument("--loss-dir", required=True, type=Path)
    parser.add_argument("--pseudogenized-suffix", default="_decayed_genes.txt")
    parser.add_argument("--deleted-suffix", default="_deleted_genes.txt")
    parser.add_argument(
        "--shared-loss-gene-list",
        type=Path,
        help=(
            "Optional exact gene list called lost in every selected sample. When supplied, "
            "write per-sample non-shared positive-call counts/rates and require every listed "
            "gene to have a positive call in every selected sample."
        ),
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    if args.expected_samples is not None and args.expected_samples < 1:
        raise SystemExit("ERROR: --expected-samples must be positive")

    try:
        reference_ids = read_reference_ids(args.reference_coords)
        reference_set = set(reference_ids)
        shared_loss_ids = (
            set(read_id_lines(args.shared_loss_gene_list, allow_empty=False))
            if args.shared_loss_gene_list else set()
        )
        outside_shared = shared_loss_ids.difference(reference_set)
        if outside_shared:
            example = ", ".join(sorted(outside_shared)[:5])
            raise ValueError(
                f"{args.shared_loss_gene_list}: {len(outside_shared)} shared-loss IDs are outside "
                f"the reference universe ({example})"
            )
        if len(shared_loss_ids) >= len(reference_ids):
            raise ValueError("Shared-loss exclusion leaves no non-shared reference genes")
        samples = read_samples(args.sample_manifest, args.include_column, args.expected_samples)
        calls: dict[str, tuple[set[str], set[str]]] = {}
        input_rows: list[dict[str, object]] = []
        for row in samples:
            sample = row["target_haplotype"]
            stem = row.get("legacy_loss_fasta_stem", "").strip() or sample
            pseudogenized_path = args.loss_dir / f"{stem}{args.pseudogenized_suffix}"
            deleted_path = args.loss_dir / f"{stem}{args.deleted_suffix}"
            pseudogenized = set(read_id_lines(pseudogenized_path, allow_empty=True))
            deleted = set(read_id_lines(deleted_path, allow_empty=True))
            overlap = pseudogenized.intersection(deleted)
            if overlap:
                raise ValueError(f"{sample}: {len(overlap)} genes occur in both positive call classes")
            outside = pseudogenized.union(deleted).difference(reference_set)
            if outside:
                example = ", ".join(sorted(outside)[:5])
                raise ValueError(f"{sample}: {len(outside)} positive calls are outside the reference universe ({example})")
            calls[sample] = (pseudogenized, deleted)
            missing_shared_positive = shared_loss_ids.difference(pseudogenized.union(deleted))
            if missing_shared_positive:
                example = ", ".join(sorted(missing_shared_positive)[:5])
                raise ValueError(
                    f"{sample}: {len(missing_shared_positive)} IDs from --shared-loss-gene-list "
                    f"lack a positive call ({example})"
                )
            for call_class, path, values in (
                ("pseudogenized", pseudogenized_path, pseudogenized),
                ("deleted", deleted_path, deleted),
            ):
                input_rows.append({
                    "target_haplotype": sample,
                    "classification": call_class,
                    "source_basename": path.name,
                    "source_sha256": sha256(path),
                    "unique_gene_ids": len(values),
                })
    except (OSError, ValueError, csv.Error) as exc:
        raise SystemExit(f"ERROR: {exc}")

    expected_rows = len(reference_ids) * len(samples)
    positive_rows = sum(len(pseudogenized) + len(deleted) for pseudogenized, deleted in calls.values())
    not_called_rows = expected_rows - positive_rows
    if not_called_rows < 0:
        raise SystemExit("ERROR: positive call count exceeds the complete matrix size")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    matrix_path = args.output_dir / "operational_loss_call_matrix.tsv.gz"
    temporary_path = args.output_dir / ".operational_loss_call_matrix.tsv.gz.tmp"
    summary_rows: list[dict[str, object]] = []
    nonshared_summary_rows: list[dict[str, object]] = []
    written = 0
    try:
        with gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["reference_gene_id", "target_haplotype", "species", "ploidy", "classification"])
            for row in samples:
                sample = row["target_haplotype"]
                pseudogenized, deleted = calls[sample]
                for gene in reference_ids:
                    classification = (
                        "pseudogenized" if gene in pseudogenized
                        else "deleted" if gene in deleted
                        else "not_called_loss"
                    )
                    writer.writerow([gene, sample, row["species"], row["ploidy"], classification])
                    written += 1
                summary_rows.append({
                    "target_haplotype": sample,
                    "species": row["species"],
                    "ploidy": row["ploidy"],
                    "reference_gene_count": len(reference_ids),
                    "pseudogenized_count": len(pseudogenized),
                    "deleted_count": len(deleted),
                    "positive_loss_call_count": len(pseudogenized) + len(deleted),
                    "not_called_loss_count": len(reference_ids) - len(pseudogenized) - len(deleted),
                })
                if shared_loss_ids:
                    nonshared_pseudogenized = len(pseudogenized.difference(shared_loss_ids))
                    nonshared_deleted = len(deleted.difference(shared_loss_ids))
                    nonshared_reference_count = len(reference_ids) - len(shared_loss_ids)
                    nonshared_positive = nonshared_pseudogenized + nonshared_deleted
                    nonshared_summary_rows.append({
                        "target_haplotype": sample,
                        "species": row["species"],
                        "ploidy": row["ploidy"],
                        "nonshared_reference_gene_count": nonshared_reference_count,
                        "excluded_shared_loss_gene_count": len(shared_loss_ids),
                        "nonshared_pseudogenized_count": nonshared_pseudogenized,
                        "nonshared_deleted_count": nonshared_deleted,
                        "nonshared_positive_loss_call_count": nonshared_positive,
                        "nonshared_not_called_loss_count": nonshared_reference_count - nonshared_positive,
                        "nonshared_positive_loss_call_rate": format(
                            nonshared_positive / nonshared_reference_count, ".12g"
                        ),
                    })
        if written != expected_rows:
            raise RuntimeError(f"wrote {written} rows; expected {expected_rows}")
        temporary_path.replace(matrix_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    with (args.output_dir / "sample_call_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    if nonshared_summary_rows:
        with (args.output_dir / "nonshared_sample_call_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(nonshared_summary_rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(nonshared_summary_rows)
    with (args.output_dir / "positive_call_input_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(input_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(input_rows)
    qc_rows = [
        ("reference_gene_count", len(reference_ids)),
        ("sample_count", len(samples)),
        ("expected_matrix_rows", expected_rows),
        ("written_matrix_rows", written),
        ("positive_loss_call_rows", positive_rows),
        ("not_called_loss_rows", not_called_rows),
    ]
    if shared_loss_ids:
        nonshared_reference_count = len(reference_ids) - len(shared_loss_ids)
        nonshared_positive_rows = sum(
            len(pseudogenized.difference(shared_loss_ids)) + len(deleted.difference(shared_loss_ids))
            for pseudogenized, deleted in calls.values()
        )
        qc_rows.extend([
            ("shared_loss_exclusion_gene_count", len(shared_loss_ids)),
            ("nonshared_reference_gene_count", nonshared_reference_count),
            ("nonshared_expected_matrix_rows", nonshared_reference_count * len(samples)),
            ("nonshared_positive_loss_call_rows", nonshared_positive_rows),
            (
                "nonshared_not_called_loss_rows",
                nonshared_reference_count * len(samples) - nonshared_positive_rows,
            ),
        ])
    with (args.output_dir / "matrix_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "value"])
        writer.writerows(qc_rows)
    metadata = [
        ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
        ("reference_coords_basename", args.reference_coords.name),
        ("reference_coords_sha256", sha256(args.reference_coords)),
        ("sample_manifest_basename", args.sample_manifest.name),
        ("sample_manifest_sha256", sha256(args.sample_manifest)),
        ("matrix_sha256", sha256(matrix_path)),
        ("absence_classification", "not_called_loss"),
        ("absence_semantics", "no positive loss call; never evidence of biological retention"),
        (
            "shared_loss_gene_list_basename",
            args.shared_loss_gene_list.name if args.shared_loss_gene_list else "",
        ),
        (
            "shared_loss_gene_list_sha256",
            sha256(args.shared_loss_gene_list) if args.shared_loss_gene_list else "",
        ),
        ("shared_loss_gene_count", len(shared_loss_ids)),
        (
            "nonshared_rate_semantics",
            "positive calls after shared-loss exclusion / non-shared reference genes; not_called_loss is not proven retention",
        ),
    ]
    with (args.output_dir / "run_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerows(metadata)
    print(f"Wrote {written} operational call rows to {matrix_path}")


if __name__ == "__main__":
    main()
