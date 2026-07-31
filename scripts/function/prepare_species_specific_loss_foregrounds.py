#!/usr/bin/env python3
"""Prepare pure single-terminal article-method loss foregrounds.

The input tree-pattern table was built from the article-comparable
``decayed + deleted`` matrix.  This adapter keeps only genes whose exact
topology placement is one single terminal-branch loss.  Recurrent independent
losses, internal-branch losses, partial/homeolog-specific states, and genes
with missing lineage states are excluded.
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
from collections import defaultdict
from pathlib import Path


class ForegroundError(ValueError):
    """Raised when frozen topology results do not support the declared set."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-tree-patterns", required=True, type=Path)
    parser.add_argument("--tree-loss-events", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-lineages", type=int, default=13)
    parser.add_argument("--expected-genes", type=int, default=1167)
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


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ForegroundError(f"missing or empty input: {path.name}")
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    if not fields or not rows or len(fields) != len(set(fields)):
        raise ForegroundError(f"invalid TSV input: {path.name}")
    return rows, fields


def require(fields: list[str], needed: set[str], label: str) -> None:
    missing = sorted(needed - set(fields))
    if missing:
        raise ForegroundError(f"{label} missing columns: {', '.join(missing)}")


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
    *,
    compressed: bool = False,
) -> None:
    opener = gzip.open if compressed else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare(
    patterns_path: Path,
    events_path: Path,
    *,
    expected_lineages: int,
    expected_genes: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    patterns, pattern_fields = read_tsv(patterns_path)
    require(
        pattern_fields,
        {
            "reference_gene_id",
            "tree_pattern",
            "tree_placement_exact",
            "complete_loss_lineage_count",
            "partial_loss_lineage_count",
            "unknown_lineage_count",
            "minimum_loss_event_count",
        },
        patterns_path.name,
    )
    selected: set[str] = set()
    for line_number, row in enumerate(patterns, 2):
        if row["tree_pattern"] != "single_terminal_branch_loss":
            continue
        try:
            exact = row["tree_placement_exact"].lower() == "true"
            complete = int(row["complete_loss_lineage_count"])
            partial = int(row["partial_loss_lineage_count"])
            unknown = int(row["unknown_lineage_count"])
            events = int(row["minimum_loss_event_count"])
        except ValueError as error:
            raise ForegroundError(
                f"{patterns_path.name}:{line_number}: invalid pattern counts"
            ) from error
        if not exact or complete != 1 or events != 1:
            raise ForegroundError(
                f"{patterns_path.name}:{line_number}: single-terminal definition mismatch"
            )
        if partial or unknown:
            continue
        gene = row["reference_gene_id"].strip()
        if not gene or gene in selected:
            raise ForegroundError(
                f"{patterns_path.name}:{line_number}: empty or duplicate selected gene"
            )
        selected.add(gene)
    if len(selected) != expected_genes:
        raise ForegroundError(
            f"observed {len(selected)} single-terminal genes; expected {expected_genes}"
        )

    events, event_fields = read_tsv(events_path)
    require(
        event_fields,
        {
            "reference_gene_id",
            "branch_id",
            "branch_type",
            "descendant_lineage_count",
            "descendant_lineages",
            "gene_minimum_loss_event_count",
        },
        events_path.name,
    )
    by_branch: dict[str, set[str]] = defaultdict(set)
    branch_species: dict[str, str] = {}
    seen: set[str] = set()
    for line_number, row in enumerate(events, 2):
        gene = row["reference_gene_id"].strip()
        if gene not in selected:
            continue
        branch = row["branch_id"].strip()
        species = row["descendant_lineages"].strip()
        try:
            descendant_count = int(row["descendant_lineage_count"])
            minimum_events = int(row["gene_minimum_loss_event_count"])
        except ValueError as error:
            raise ForegroundError(
                f"{events_path.name}:{line_number}: invalid event counts"
            ) from error
        if (
            not branch.startswith("terminal__")
            or row["branch_type"] != "terminal"
            or descendant_count != 1
            or ";" in species
            or minimum_events != 1
        ):
            raise ForegroundError(
                f"{events_path.name}:{line_number}: selected gene is not one terminal event"
            )
        if gene in seen:
            raise ForegroundError(
                f"{events_path.name}:{line_number}: selected gene has duplicate events"
            )
        seen.add(gene)
        if branch in branch_species and branch_species[branch] != species:
            raise ForegroundError(f"{events_path.name}: branch lineage changed")
        branch_species[branch] = species
        by_branch[branch].add(gene)
    if seen != selected:
        raise ForegroundError("selected pattern genes and terminal event genes differ")
    if len(by_branch) != expected_lineages or len(set(branch_species.values())) != expected_lineages:
        raise ForegroundError(
            f"observed {len(by_branch)} terminal lineage foregrounds; expected {expected_lineages}"
        )

    membership_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    for branch in sorted(by_branch, key=lambda value: branch_species[value]):
        foreground_id = f"species_specific__{branch}"
        genes = sorted(by_branch[branch])
        membership_rows.extend(
            {"foreground_id": foreground_id, "reference_gene_id": gene}
            for gene in genes
        )
        metadata_rows.append(
            {
                "foreground_id": foreground_id,
                "analysis_scope": "single_terminal_branch_species_specific_loss",
                "background_scope": "all_13_lineages_resolved_article_method",
                "branch_id": branch,
                "descendant_lineage_count": 1,
                "descendant_lineages": branch_species[branch],
                "foreground_gene_count": len(genes),
                "definition": (
                    "article-method decayed+deleted; exact one terminal event; "
                    "all other lineages resolved and retained"
                ),
            }
        )
    if sum(int(row["foreground_gene_count"]) for row in metadata_rows) != expected_genes:
        raise ForegroundError("foreground counts do not close to selected genes")
    summary = {
        "schema_version": "1.0",
        "status": "PASS_ARTICLE_METHOD_SINGLE_TERMINAL_FOREGROUNDS",
        "loss_classification": "article-method decayed + deleted",
        "species_specific_definition": (
            "one exact terminal-branch event only; recurrent, internal, partial, "
            "and unknown patterns excluded"
        ),
        "lineage_foregrounds": len(metadata_rows),
        "species_specific_loss_genes": len(selected),
    }
    return membership_rows, metadata_rows, summary


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise ForegroundError(f"output directory already exists: {args.output_dir}")
    membership, metadata, summary = prepare(
        args.gene_tree_patterns,
        args.tree_loss_events,
        expected_lineages=args.expected_lineages,
        expected_genes=args.expected_genes,
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
        metadata_path = staging / "foreground_metadata.tsv"
        write_tsv(
            genes_path,
            ["foreground_id", "reference_gene_id"],
            membership,
            compressed=True,
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
                "definition",
            ],
            metadata,
        )
        summary["inputs"] = [
            {
                "role": "gene_tree_patterns",
                "basename": args.gene_tree_patterns.name,
                "sha256": sha256(args.gene_tree_patterns),
            },
            {
                "role": "tree_loss_events",
                "basename": args.tree_loss_events.name,
                "sha256": sha256(args.tree_loss_events),
            },
        ]
        summary["outputs"] = [
            {
                "basename": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in (genes_path, metadata_path)
        ]
        (staging / "run_manifest.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
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
