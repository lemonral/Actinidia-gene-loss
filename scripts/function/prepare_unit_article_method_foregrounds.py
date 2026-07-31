#!/usr/bin/env python3
"""Prepare 23 independent assembly-unit loss foregrounds.

Each foreground uses the article-method positive rule (``decayed + deleted``)
for one assembly unit.  Units, haplotypes, and subgenomes are never aggregated.
The matching enrichment background contains the resolved article-method rows
(``retained + decayed + deleted``) for that same unit; ``not_called_loss`` is
excluded from both the foreground and its background.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


POSITIVE = {"decayed", "deleted"}
RESOLVED = {"retained", "decayed", "deleted"}
ALL_CLASSES = {*RESOLVED, "not_called_loss"}


class ForegroundError(ValueError):
    """Raised when unit-level article-method inputs are inconsistent."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-matrix", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-matrix-sha256", default="")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def read_fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].split(None, 1)[0]
            if not identifier or identifier in seen:
                raise ForegroundError(
                    f"{path.name}:{line_number}: empty or duplicate FASTA ID"
                )
            seen.add(identifier)
            identifiers.append(identifier)
    if not identifiers:
        raise ForegroundError(f"{path.name}: no FASTA records")
    return identifiers


def read_metadata(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "include",
        }
        if not required.issubset(fields):
            raise ForegroundError(f"{path.name}: missing unit metadata columns")
        rows = [
            dict(row)
            for row in reader
            if row["include"].strip().lower() == "true"
        ]
    units: list[str] = []
    metadata: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, 2):
        unit = row["assembly_unit_id"].strip()
        species = row["biological_species"].strip()
        if not unit or not species or unit in metadata:
            raise ForegroundError(
                f"{path.name}:{line_number}: empty or duplicate included unit"
            )
        units.append(unit)
        metadata[unit] = row
    return units, metadata


def write_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[Mapping[str, object]],
    *,
    gz: bool = False,
) -> None:
    opener = gzip.open if gz else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise ForegroundError(f"output directory already exists: {args.output_dir}")
    for path in (args.unit_matrix, args.unit_metadata, args.reference_protein):
        if not path.is_file() or path.stat().st_size == 0:
            raise ForegroundError(f"missing or empty input: {path}")
    matrix_hash = sha256(args.unit_matrix)
    if (
        args.expected_matrix_sha256
        and matrix_hash != args.expected_matrix_sha256.lower()
    ):
        raise ForegroundError("unit matrix SHA-256 mismatch")

    reference_order = read_fasta_ids(args.reference_protein)
    if len(reference_order) != args.expected_reference_genes:
        raise ForegroundError(
            f"observed {len(reference_order)} reference genes; "
            f"expected {args.expected_reference_genes}"
        )
    reference = set(reference_order)
    units, metadata = read_metadata(args.unit_metadata)
    if len(units) != args.expected_units:
        raise ForegroundError(
            f"observed {len(units)} included units; expected {args.expected_units}"
        )

    foregrounds: dict[str, set[str]] = defaultdict(set)
    backgrounds: dict[str, set[str]] = defaultdict(set)
    class_counts: dict[str, Counter[str]] = defaultdict(Counter)
    seen: set[tuple[str, str]] = set()
    with open_text(args.unit_matrix) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "reference_gene_id",
            "assembly_unit_id",
            "manuscript_classification",
            "manuscript_positive_loss",
        }
        if not required.issubset(fields):
            raise ForegroundError(
                f"{args.unit_matrix.name}: missing required matrix columns"
            )
        for line_number, row in enumerate(reader, 2):
            gene = row["reference_gene_id"].strip()
            unit = row["assembly_unit_id"].strip()
            state = row["manuscript_classification"].strip()
            positive = row["manuscript_positive_loss"].strip().lower()
            pair = (unit, gene)
            if (
                unit not in metadata
                or gene not in reference
                or state not in ALL_CLASSES
                or positive not in {"true", "false"}
                or (positive == "true") != (state in POSITIVE)
                or pair in seen
            ):
                raise ForegroundError(
                    f"{args.unit_matrix.name}:{line_number}: invalid unit-gene row"
                )
            seen.add(pair)
            class_counts[unit][state] += 1
            if state in RESOLVED:
                backgrounds[unit].add(gene)
            if state in POSITIVE:
                foregrounds[unit].add(gene)

    expected_rows = args.expected_units * args.expected_reference_genes
    if len(seen) != expected_rows:
        raise ForegroundError(
            f"observed {len(seen)} unit-gene rows; expected {expected_rows}"
        )
    if any(
        len({gene for candidate, gene in seen if candidate == unit})
        != args.expected_reference_genes
        for unit in units
    ):
        raise ForegroundError("unit matrix is not a complete rectangular grid")
    if any(not foregrounds[unit] for unit in units):
        raise ForegroundError("at least one assembly-unit foreground is empty")

    foreground_rows: list[dict[str, object]] = []
    background_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for unit in units:
        foreground_id = f"unit__{unit}"
        background_scope = f"resolved__{unit}"
        for gene in sorted(foregrounds[unit]):
            foreground_rows.append(
                {
                    "foreground_id": foreground_id,
                    "reference_gene_id": gene,
                }
            )
        for gene in sorted(backgrounds[unit]):
            background_rows.append(
                {
                    "background_scope": background_scope,
                    "reference_gene_id": gene,
                }
            )
        species = metadata[unit]["biological_species"].strip()
        suffix = metadata[unit]["haplotype_or_subgenome"].strip()
        metadata_rows.append(
            {
                "foreground_id": foreground_id,
                "analysis_scope": "assembly_unit_article_method_loss",
                "background_scope": background_scope,
                "branch_id": unit,
                "descendant_lineage_count": 1,
                "descendant_lineages": species,
                "foreground_gene_count": len(foregrounds[unit]),
                "assembly_unit_id": unit,
                "biological_species": species,
                "haplotype_or_subgenome": suffix,
            }
        )
        summary_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": species,
                "haplotype_or_subgenome": suffix,
                "retained": class_counts[unit]["retained"],
                "decayed": class_counts[unit]["decayed"],
                "deleted": class_counts[unit]["deleted"],
                "not_called_loss": class_counts[unit]["not_called_loss"],
                "positive_loss": len(foregrounds[unit]),
                "resolved_background": len(backgrounds[unit]),
                "loss_rate": (
                    f"{len(foregrounds[unit]) / len(backgrounds[unit]):.12f}"
                ),
            }
        )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        genes_path = staging / "foreground_gene_ids.tsv.gz"
        backgrounds_path = staging / "foreground_background_gene_ids.tsv.gz"
        metadata_path = staging / "foreground_metadata.tsv"
        summary_path = staging / "unit_loss_summary.tsv"
        write_tsv(
            genes_path,
            ["foreground_id", "reference_gene_id"],
            foreground_rows,
            gz=True,
        )
        write_tsv(
            backgrounds_path,
            ["background_scope", "reference_gene_id"],
            background_rows,
            gz=True,
        )
        write_tsv(
            metadata_path,
            [
                "foreground_id",
                "analysis_scope",
                "background_scope",
                "branch_id",
                "descendant_lineage_count",
                "descendant_lineages",
                "foreground_gene_count",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
            ],
            metadata_rows,
        )
        write_tsv(
            summary_path,
            list(summary_rows[0]),
            summary_rows,
        )
        output_paths = [
            genes_path,
            backgrounds_path,
            metadata_path,
            summary_path,
        ]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_UNIT_ARTICLE_METHOD_FOREGROUNDS",
            "created_at_utc": datetime.now(timezone.utc).isoformat(
                timespec="seconds"
            ),
            "definitions": {
                "foreground": (
                    "decayed + deleted within one assembly unit"
                ),
                "background": (
                    "retained + decayed + deleted within the same assembly "
                    "unit; not_called_loss excluded"
                ),
                "aggregation": (
                    "none; all haplotypes and subgenomes remain independent"
                ),
            },
            "counts": {
                "assembly_units": len(units),
                "reference_genes": len(reference),
                "matrix_rows": len(seen),
                "foreground_memberships": sum(
                    len(foregrounds[unit]) for unit in units
                ),
                "resolved_background_memberships": sum(
                    len(backgrounds[unit]) for unit in units
                ),
            },
            "inputs": [
                {
                    "role": role,
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for role, path in (
                    ("unit_matrix", args.unit_matrix),
                    ("unit_metadata", args.unit_metadata),
                    ("reference_protein", args.reference_protein),
                )
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in output_paths
            ],
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, args.output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (ForegroundError, OSError, csv.Error, UnicodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
