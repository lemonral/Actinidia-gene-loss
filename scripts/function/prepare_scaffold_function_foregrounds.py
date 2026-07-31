#!/usr/bin/env python3
"""Prepare GO/KEGG foregrounds for the 23-unit topology scaffold."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


class ForegroundError(ValueError):
    """Raised when scaffold events cannot form exact enrichment inputs."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaffold-events", required=True, type=Path)
    parser.add_argument("--gene-patterns", required=True, type=Path)
    parser.add_argument("--scaffold-nodes", required=True, type=Path)
    parser.add_argument("--scaffold-manifest", required=True, type=Path)
    parser.add_argument("--expected-nodes", type=int, default=39)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-resolved-genes", type=int, default=33998)
    parser.add_argument("--expected-event-rows", type=int, default=56602)
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
    return path.open(encoding="utf-8-sig", newline="")


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = reader.fieldnames or []
    if not rows or len(fields) != len(set(fields)):
        raise ForegroundError(f"{path.name}: invalid or empty table")
    return rows, fields


def require(fields: Iterable[str], needed: set[str], label: str) -> None:
    missing = sorted(needed - set(fields))
    if missing:
        raise ForegroundError(f"{label}: missing {', '.join(missing)}")


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
        writer.writerows(rows)


def output_record(path: Path) -> dict[str, object]:
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise ForegroundError(f"output directory exists: {args.output_dir}")
    inputs = {
        "scaffold_events": args.scaffold_events,
        "gene_patterns": args.gene_patterns,
        "scaffold_nodes": args.scaffold_nodes,
        "scaffold_manifest": args.scaffold_manifest,
    }
    for role, path in inputs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise ForegroundError(f"missing or empty {role}: {path}")

    scaffold_manifest = json.loads(
        args.scaffold_manifest.read_text(encoding="utf-8")
    )
    if (
        scaffold_manifest.get("status")
        != "PASS_UNIT_RESOLVED_ARTICLE_LOSS_SCAFFOLD"
    ):
        raise ForegroundError("scaffold manifest is not PASS")

    nodes, node_fields = read_tsv(args.scaffold_nodes)
    require(
        node_fields,
        {
            "node_id",
            "parent_node_id",
            "node_type",
            "node_name",
            "descendant_unit_count",
            "descendant_units",
            "minimum_leaf_plot_order",
        },
        args.scaffold_nodes.name,
    )
    if len(nodes) != args.expected_nodes:
        raise ForegroundError(
            f"expected {args.expected_nodes} nodes, observed {len(nodes)}"
        )
    node_by_id = {row["node_id"]: row for row in nodes}
    if len(node_by_id) != len(nodes):
        raise ForegroundError("duplicate scaffold node")

    patterns, pattern_fields = read_tsv(args.gene_patterns)
    require(
        pattern_fields,
        {"reference_gene_id", "scaffold_pattern"},
        args.gene_patterns.name,
    )
    if len(patterns) != args.expected_reference_genes:
        raise ForegroundError("reference-gene count mismatch")
    pattern_by_gene = {row["reference_gene_id"]: row for row in patterns}
    if len(pattern_by_gene) != len(patterns):
        raise ForegroundError("duplicate gene pattern")
    resolved = sorted(
        gene
        for gene, row in pattern_by_gene.items()
        if row["scaffold_pattern"] != "ambiguous_not_called"
    )
    if len(resolved) != args.expected_resolved_genes:
        raise ForegroundError(
            f"expected {args.expected_resolved_genes} resolved genes, "
            f"observed {len(resolved)}"
        )
    resolved_set = set(resolved)

    events, event_fields = read_tsv(args.scaffold_events)
    require(
        event_fields,
        {"reference_gene_id", "node_id"},
        args.scaffold_events.name,
    )
    if len(events) != args.expected_event_rows:
        raise ForegroundError(
            f"expected {args.expected_event_rows} event rows, "
            f"observed {len(events)}"
        )
    seen: set[tuple[str, str]] = set()
    event_counts: Counter[str] = Counter()
    membership: list[dict[str, str]] = []
    for line_number, row in enumerate(events, 2):
        gene = row["reference_gene_id"]
        node_id = row["node_id"]
        key = (node_id, gene)
        if (
            not gene
            or node_id not in node_by_id
            or gene not in resolved_set
            or key in seen
        ):
            raise ForegroundError(
                f"{args.scaffold_events.name}:{line_number}: invalid event"
            )
        seen.add(key)
        event_counts[node_id] += 1
        membership.append(
            {
                "foreground_id": f"scaffold__{node_id}",
                "reference_gene_id": gene,
            }
        )
    if set(event_counts) != set(node_by_id):
        raise ForegroundError("one or more scaffold nodes lack an event set")

    metadata: list[dict[str, object]] = []
    for node in sorted(
        nodes,
        key=lambda row: (
            int(row["minimum_leaf_plot_order"]),
            int(row["descendant_unit_count"]),
            row["node_id"],
        ),
    ):
        node_id = node["node_id"]
        metadata.append(
            {
                "foreground_id": f"scaffold__{node_id}",
                "analysis_scope": "assembly_unit_topology_scaffold_loss",
                "background_scope": (
                    "all_23_units_resolved_article_method"
                ),
                "foreground_gene_count": event_counts[node_id],
                "branch_id": node_id,
                "descendant_lineage_count": node[
                    "descendant_unit_count"
                ],
                "descendant_lineages": node["descendant_units"],
                "node_type": node["node_type"],
                "node_name": node["node_name"],
                "parent_node_id": node["parent_node_id"],
                "minimum_leaf_plot_order": node[
                    "minimum_leaf_plot_order"
                ],
            }
        )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.",
            dir=args.output_dir.parent,
        )
    )
    try:
        membership_path = temporary / "foreground_gene_ids.tsv.gz"
        metadata_path = temporary / "foreground_metadata.tsv"
        background_path = temporary / "resolved_background_gene_ids.txt"
        write_tsv(
            membership_path,
            ["foreground_id", "reference_gene_id"],
            membership,
            gz=True,
        )
        write_tsv(
            metadata_path,
            list(metadata[0]),
            metadata,
        )
        background_path.write_text(
            "".join(f"{gene}\n" for gene in resolved),
            encoding="utf-8",
        )
        outputs = [membership_path, metadata_path, background_path]
        manifest = {
            "schema_version": 1,
            "status": "PASS_UNIT_SCAFFOLD_FUNCTION_FOREGROUNDS",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "definitions": {
                "foreground": (
                    "maximal decayed-plus-deleted loss-event genes assigned "
                    "to one node of the 23-unit topology scaffold"
                ),
                "background": (
                    "reference genes with resolved article-method states in "
                    "all 23 assembly units"
                ),
                "topology_claim": (
                    "topology-only assembly-unit scaffold; not a newly "
                    "inferred 23-species phylogeny"
                ),
            },
            "counts": {
                "foregrounds": len(metadata),
                "foreground_memberships": len(membership),
                "resolved_background_genes": len(resolved),
                "reference_genes": len(patterns),
            },
            "inputs": [
                {
                    "role": role,
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for role, path in inputs.items()
            ],
            "outputs": [output_record(path) for path in outputs],
        }
        manifest_path = temporary / "run_manifest.json"
        manifest_path.write_text(
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
    except (ForegroundError, OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
