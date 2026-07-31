#!/usr/bin/env python3
"""Prepare unit-resolved article-method loss evidence and a 23-tip scaffold.

The publication loss definition is the historical ``decayed + deleted`` rule.
The refined coding-disruption fields are annotations of that fixed result and
never rewrite the historical class.  Assembly units remain separate leaves.
For multi-unit biological species, leaves are attached as a polytomy at the
matching species tip of the accepted 13-lineage backbone.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ARTICLE_CLASSES = {"retained", "decayed", "deleted", "not_called_loss"}
POSITIVE = {"decayed", "deleted"}
CONFIRMED_CAUSES = {
    "frameshift_supported": "frameshift_only",
    "stop_supported": "inframe_stop_only",
    "frameshift_and_stop_supported": "frameshift_and_inframe_stop",
}
CANDIDATE_CAUSES = {
    "n_terminal_alignment_truncation_candidate",
    "c_terminal_alignment_truncation_candidate",
    "both_terminal_alignment_truncation_candidate",
    "partial_local_alignment_other_candidate",
}


class ScaffoldError(ValueError):
    """Raised when frozen inputs cannot support the declared scaffold."""


@dataclass
class Node:
    """Minimal rooted tree node used for topology-only event placement."""

    name: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    node_id: str = ""
    node_type: str = ""

    def leaves(self) -> list["Node"]:
        if not self.children:
            return [self]
        return [leaf for child in self.children for leaf in child.leaves()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-matrix", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--tip-map", required=True, type=Path)
    parser.add_argument("--time-tree", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-lineages", type=int, default=13)
    parser.add_argument("--expected-positive-rows", type=int, default=179827)
    parser.add_argument("--expected-shared-genes", type=int, default=3616)
    parser.add_argument("--expected-decayed-frameshift-only", type=int, default=11559)
    parser.add_argument("--expected-decayed-stop-only", type=int, default=3258)
    parser.add_argument("--expected-decayed-both", type=int, default=5071)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, *, usecols: list[str] | None = None) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise ScaffoldError(f"missing or empty input: {path.name}")
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        compression="infer",
        usecols=usecols,
    )


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ScaffoldError(f"{label} missing columns: {', '.join(missing)}")


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def tokenize_newick(text: str) -> list[str]:
    tokens = re.findall(r"\(|\)|,|:|;|[^\s(),:;]+", text.strip())
    if not tokens or tokens[-1] != ";":
        raise ScaffoldError("time tree must be one semicolon-terminated Newick")
    return tokens


def parse_newick(path: Path) -> Node:
    tokens = tokenize_newick(path.read_text(encoding="utf-8"))
    index = 0

    def subtree() -> Node:
        nonlocal index
        if index >= len(tokens):
            raise ScaffoldError("unexpected end of Newick tree")
        if tokens[index] == "(":
            index += 1
            children = [subtree()]
            while index < len(tokens) and tokens[index] == ",":
                index += 1
                children.append(subtree())
            if index >= len(tokens) or tokens[index] != ")":
                raise ScaffoldError("unbalanced Newick parentheses")
            index += 1
            name = ""
            if index < len(tokens) and tokens[index] not in {":", ",", ")", ";"}:
                name = tokens[index]
                index += 1
            node = Node(name=name, children=children)
            for child in children:
                child.parent = node
        else:
            if tokens[index] in {",", ")", ":", ";"}:
                raise ScaffoldError(f"invalid Newick token: {tokens[index]!r}")
            node = Node(name=tokens[index])
            index += 1
        if index < len(tokens) and tokens[index] == ":":
            index += 1
            if index >= len(tokens):
                raise ScaffoldError("missing Newick branch length")
            try:
                float(tokens[index])
            except ValueError as error:
                raise ScaffoldError(
                    f"invalid Newick branch length: {tokens[index]!r}"
                ) from error
            index += 1
        return node

    root = subtree()
    if index != len(tokens) - 1:
        raise ScaffoldError("unexpected tokens after Newick root")
    names = [leaf.name for leaf in root.leaves()]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ScaffoldError("Newick tips must be nonempty and unique")
    return root


def walk(node: Node) -> list[Node]:
    return [node, *(item for child in node.children for item in walk(child))]


def descendant_names(node: Node) -> set[str]:
    return {leaf.name for leaf in node.leaves()}


def read_metadata(
    path: Path,
    expected_units: int,
    expected_lineages: int,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    frame = read_tsv(path)
    require_columns(
        frame,
        {
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "include",
        },
        "unit metadata",
    )
    frame = frame.loc[frame["include"].str.lower() == "true"].copy()
    if len(frame) != expected_units or frame["assembly_unit_id"].duplicated().any():
        raise ScaffoldError("unit metadata is not the exact unique unit cohort")
    species_units: dict[str, list[str]] = defaultdict(list)
    for row in frame.itertuples(index=False):
        species_units[row.biological_species].append(row.assembly_unit_id)
    if len(species_units) != expected_lineages:
        raise ScaffoldError("unit metadata lineage count changed")
    return frame, dict(species_units)


def actinidia_subtree(
    root: Node,
    tip_map_path: Path,
    species: set[str],
) -> tuple[Node, dict[str, str]]:
    frame = read_tsv(tip_map_path)
    require_columns(
        frame,
        {"tree_tip", "biological_species", "include"},
        "tip map",
    )
    frame = frame.loc[frame["include"].str.lower() == "true"].copy()
    if frame["tree_tip"].duplicated().any() or frame[
        "biological_species"
    ].duplicated().any():
        raise ScaffoldError("included tip map must be one-to-one")
    tip_to_species = dict(zip(frame["tree_tip"], frame["biological_species"]))
    if set(tip_to_species.values()) != species:
        raise ScaffoldError("tip-map lineages do not match unit metadata")
    tree_tips = descendant_names(root)
    if not set(tip_to_species).issubset(tree_tips):
        raise ScaffoldError("included tip is absent from the time tree")
    wanted = set(tip_to_species)
    candidates = [
        node for node in walk(root) if wanted.issubset(descendant_names(node))
    ]
    subtree = min(candidates, key=lambda node: len(descendant_names(node)))
    if descendant_names(subtree) != wanted:
        raise ScaffoldError("included tips do not form one exact Actinidia clade")
    subtree.parent = None
    return subtree, tip_to_species


def expand_units(
    root: Node,
    tip_to_species: dict[str, str],
    species_units: dict[str, list[str]],
) -> None:
    original_leaves = list(root.leaves())
    for leaf in original_leaves:
        species = tip_to_species[leaf.name]
        units = species_units[species]
        if len(units) == 1:
            leaf.name = units[0]
            leaf.node_type = "unit_terminal"
            continue
        leaf.name = species
        leaf.node_type = "species_unit_group"
        leaf.children = [
            Node(name=unit, parent=leaf, node_type="unit_terminal")
            for unit in units
        ]


def assign_node_ids(root: Node) -> dict[str, Node]:
    by_id: dict[str, Node] = {}
    for node in walk(root):
        descendants = sorted(leaf.name for leaf in node.leaves())
        if not node.children:
            node.node_type = "unit_terminal"
            node.node_id = f"unit__{node.name}"
        elif node.node_type == "species_unit_group":
            node.node_id = f"group__{slug(node.name)}"
        elif node is root:
            node.node_type = "backbone_internal"
            node.node_id = "internal__actinidia_all"
        else:
            node.node_type = "backbone_internal"
            token = hashlib.sha1(
                "\n".join(descendants).encode("utf-8")
            ).hexdigest()[:10]
            node.node_id = f"internal__{token}"
        if node.node_id in by_id:
            raise ScaffoldError(f"duplicate scaffold node ID: {node.node_id}")
        by_id[node.node_id] = node
    return by_id


def node_depth(node: Node) -> int:
    depth = 0
    while node.parent is not None:
        depth += 1
        node = node.parent
    return depth


def newick_record(node: Node) -> str:
    if node.children:
        return f"({','.join(newick_record(child) for child in node.children)}){node.node_id}"
    return node.node_id


def maximal_loss_nodes(root: Node, positive_units: set[str]) -> list[Node]:
    events: list[Node] = []

    def visit(node: Node, parent_all_positive: bool) -> None:
        descendants = {leaf.name for leaf in node.leaves()}
        all_positive = bool(descendants) and descendants.issubset(positive_units)
        if all_positive and not parent_all_positive:
            events.append(node)
            return
        for child in node.children:
            visit(child, all_positive)

    visit(root, False)
    return events


def write_frame(path: Path, frame: pd.DataFrame, *, gzip_output: bool = False) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "mtime": 0} if gzip_output else None,
    )


def run(args: argparse.Namespace) -> None:
    inputs = [
        args.manuscript_matrix,
        args.unit_metadata,
        args.tip_map,
        args.time_tree,
    ]
    missing = [path.name for path in inputs if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ScaffoldError(f"missing or empty inputs: {missing}")
    if args.output_dir.exists():
        raise ScaffoldError(f"output directory already exists: {args.output_dir}")

    metadata, species_units = read_metadata(
        args.unit_metadata,
        args.expected_units,
        args.expected_lineages,
    )
    unit_order = list(metadata["assembly_unit_id"])
    unit_set = set(unit_order)
    root_all = parse_newick(args.time_tree)
    root, tip_to_species = actinidia_subtree(
        root_all,
        args.tip_map,
        set(species_units),
    )
    expand_units(root, tip_to_species, species_units)
    nodes = assign_node_ids(root)
    if len(root.leaves()) != args.expected_units or {
        leaf.name for leaf in root.leaves()
    } != unit_set:
        raise ScaffoldError("expanded scaffold does not close to the unit cohort")

    needed = [
        "reference_gene_id",
        "assembly_unit_id",
        "manuscript_classification",
        "manuscript_positive_loss",
        "uniform_classification",
        "uniform_evidence_reason",
        "refined_decayed_cause",
        "refined_cause_evidence_level",
        "query_coverage",
        "exact_alignment_identity",
        "alignment_score",
        "frameshift_events",
        "inframe_stop_codons",
    ]
    matrix = read_tsv(args.manuscript_matrix, usecols=needed)
    expected_rows = args.expected_units * args.expected_reference_genes
    if len(matrix) != expected_rows:
        raise ScaffoldError(
            f"matrix rows={len(matrix):,}; expected {expected_rows:,}"
        )
    if matrix.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
        raise ScaffoldError("matrix contains duplicate gene-unit rows")
    if set(matrix["assembly_unit_id"]) != unit_set:
        raise ScaffoldError("matrix unit universe differs from metadata")
    if set(matrix["manuscript_classification"]) - ARTICLE_CLASSES:
        raise ScaffoldError("matrix contains an unsupported historical class")
    per_gene = matrix.groupby("reference_gene_id", sort=False).size()
    if len(per_gene) != args.expected_reference_genes or set(per_gene) != {
        args.expected_units
    }:
        raise ScaffoldError("matrix does not contain one complete unit grid per gene")
    matrix["positive"] = matrix["manuscript_classification"].isin(POSITIVE)
    declared_positive = matrix["manuscript_positive_loss"].str.lower() == "true"
    if not matrix["positive"].equals(declared_positive):
        raise ScaffoldError("historical class and positive flag disagree")
    positive_rows = int(matrix["positive"].sum())
    if positive_rows != args.expected_positive_rows:
        raise ScaffoldError(
            f"positive rows={positive_rows:,}; expected {args.expected_positive_rows:,}"
        )

    positive_by_gene = matrix.groupby("reference_gene_id", sort=False)[
        "positive"
    ].sum()
    shared_genes = set(
        positive_by_gene.index[positive_by_gene == args.expected_units]
    )
    if len(shared_genes) != args.expected_shared_genes:
        raise ScaffoldError(
            f"shared genes={len(shared_genes):,}; expected {args.expected_shared_genes:,}"
        )
    matrix["shared_positive"] = (
        matrix["positive"] & matrix["reference_gene_id"].isin(shared_genes)
    )

    unit_rows: list[dict[str, object]] = []
    shared_rows: list[dict[str, object]] = []
    mechanism_rows: list[dict[str, object]] = []
    cause_rows: list[dict[str, object]] = []
    metadata_lookup = metadata.set_index("assembly_unit_id").to_dict("index")
    for unit in unit_order:
        subset = matrix.loc[matrix["assembly_unit_id"] == unit]
        counts = Counter(subset["manuscript_classification"])
        if sum(counts.values()) != args.expected_reference_genes:
            raise ScaffoldError(f"{unit}: unit rows do not close")
        positive = counts["decayed"] + counts["deleted"]
        unit_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": metadata_lookup[unit]["biological_species"],
                "haplotype_or_subgenome": metadata_lookup[unit][
                    "haplotype_or_subgenome"
                ],
                "retained": counts["retained"],
                "decayed": counts["decayed"],
                "deleted": counts["deleted"],
                "not_called_loss": counts["not_called_loss"],
                "positive_loss": positive,
                "resolved_denominator": counts["retained"] + positive,
                "positive_loss_rate": positive / (counts["retained"] + positive),
            }
        )

        positive_subset = subset.loc[subset["positive"]].copy()
        shared_class = Counter(
            zip(
                positive_subset["shared_positive"],
                positive_subset["manuscript_classification"],
            )
        )
        shared_row = {
            "assembly_unit_id": unit,
            "biological_species": metadata_lookup[unit]["biological_species"],
            "haplotype_or_subgenome": metadata_lookup[unit][
                "haplotype_or_subgenome"
            ],
            "shared_decayed": shared_class[(True, "decayed")],
            "shared_deleted": shared_class[(True, "deleted")],
            "nonshared_decayed": shared_class[(False, "decayed")],
            "nonshared_deleted": shared_class[(False, "deleted")],
        }
        shared_row["shared_positive_loss"] = (
            shared_row["shared_decayed"] + shared_row["shared_deleted"]
        )
        shared_row["nonshared_positive_loss"] = (
            shared_row["nonshared_decayed"]
            + shared_row["nonshared_deleted"]
        )
        if (
            shared_row["shared_positive_loss"]
            + shared_row["nonshared_positive_loss"]
            != positive
        ):
            raise ScaffoldError(f"{unit}: shared/non-shared counts do not close")
        shared_rows.append(shared_row)

        decayed_subset = subset.loc[
            subset["manuscript_classification"] == "decayed"
        ]
        causes = Counter(decayed_subset["refined_decayed_cause"])
        confirmed = {
            label: causes[cause] for cause, label in CONFIRMED_CAUSES.items()
        }
        confirmed_total = sum(confirmed.values())
        candidate_total = sum(causes[cause] for cause in CANDIDATE_CAUSES)
        unresolved = counts["decayed"] - confirmed_total - candidate_total
        if unresolved < 0:
            raise ScaffoldError(f"{unit}: decayed mechanism counts exceed decayed")
        mechanism_rows.append(
            {
                "assembly_unit_id": unit,
                "biological_species": metadata_lookup[unit]["biological_species"],
                "haplotype_or_subgenome": metadata_lookup[unit][
                    "haplotype_or_subgenome"
                ],
                "decayed": counts["decayed"],
                "frameshift_only": confirmed["frameshift_only"],
                "inframe_stop_only": confirmed["inframe_stop_only"],
                "frameshift_and_inframe_stop": confirmed[
                    "frameshift_and_inframe_stop"
                ],
                "confirmed_type_total": confirmed_total,
                "confirmed_type_fraction_of_decayed": confirmed_total
                / counts["decayed"],
                "candidate_type_total": candidate_total,
                "unresolved_decayed": unresolved,
            }
        )
        for cause, count in sorted(causes.items()):
            if cause in CONFIRMED_CAUSES:
                evidence_tier = "confirmed_coding_disruption"
                display_class = CONFIRMED_CAUSES[cause]
            elif cause in CANDIDATE_CAUSES:
                evidence_tier = "candidate"
                display_class = cause
            else:
                evidence_tier = "unresolved"
                display_class = "unresolved_decayed"
            cause_rows.append(
                {
                    "assembly_unit_id": unit,
                    "refined_decayed_cause": cause,
                    "evidence_tier": evidence_tier,
                    "display_class": display_class,
                    "decayed_unit_gene_rows": count,
                }
            )

    unit_frame = pd.DataFrame(unit_rows)
    shared_frame = pd.DataFrame(shared_rows)
    mechanism_frame = pd.DataFrame(mechanism_rows)
    cause_frame = pd.DataFrame(cause_rows)
    confirmed_totals = mechanism_frame[
        [
            "frameshift_only",
            "inframe_stop_only",
            "frameshift_and_inframe_stop",
        ]
    ].sum()
    expected_confirmed = {
        "frameshift_only": args.expected_decayed_frameshift_only,
        "inframe_stop_only": args.expected_decayed_stop_only,
        "frameshift_and_inframe_stop": args.expected_decayed_both,
    }
    if confirmed_totals.to_dict() != expected_confirmed:
        raise ScaffoldError(
            f"decayed confirmed-type totals changed: {confirmed_totals.to_dict()}"
        )

    mechanism_display: list[str] = []
    for row in matrix.itertuples(index=False):
        if not row.positive:
            mechanism_display.append("")
        elif row.manuscript_classification == "deleted":
            mechanism_display.append("not_applicable_deleted")
        elif row.refined_decayed_cause in CONFIRMED_CAUSES:
            mechanism_display.append(CONFIRMED_CAUSES[row.refined_decayed_cause])
        elif row.refined_decayed_cause in CANDIDATE_CAUSES:
            mechanism_display.append(row.refined_decayed_cause)
        else:
            mechanism_display.append("unresolved_decayed")
    matrix["mechanism_display_class"] = mechanism_display
    crosswalk_columns = [
        "reference_gene_id",
        "assembly_unit_id",
        "manuscript_classification",
        "uniform_classification",
        "uniform_evidence_reason",
        "refined_decayed_cause",
        "refined_cause_evidence_level",
        "mechanism_display_class",
        "query_coverage",
        "exact_alignment_identity",
        "alignment_score",
        "frameshift_events",
        "inframe_stop_codons",
    ]
    crosswalk = matrix.loc[matrix["positive"], crosswalk_columns].copy()
    if len(crosswalk) != args.expected_positive_rows:
        raise ScaffoldError("positive mechanism crosswalk does not close")

    state_grid = matrix.pivot(
        index="reference_gene_id",
        columns="assembly_unit_id",
        values="manuscript_classification",
    )[unit_order]
    strict_gene_units: dict[str, set[str]] = defaultdict(set)
    confirmed_matrix = matrix.loc[
        (matrix["manuscript_classification"] == "decayed")
        & matrix["refined_decayed_cause"].isin(CONFIRMED_CAUSES)
    ]
    for gene, unit in zip(
        confirmed_matrix["reference_gene_id"],
        confirmed_matrix["assembly_unit_id"],
    ):
        strict_gene_units[gene].add(unit)

    node_descendants = {
        node.node_id: tuple(leaf.name for leaf in node.leaves())
        for node in walk(root)
    }
    event_rows: list[dict[str, object]] = []
    pattern_rows: list[dict[str, object]] = []
    ambiguous_rows: list[dict[str, object]] = []
    pattern_counts: Counter[str] = Counter()
    for gene, row in state_grid.iterrows():
        states = row.to_dict()
        unknown_units = {
            unit for unit, state in states.items() if state == "not_called_loss"
        }
        positive_units = {unit for unit, state in states.items() if state in POSITIVE}
        if unknown_units:
            pattern = "ambiguous_not_called"
            events: list[Node] = []
            ambiguous_rows.append(
                {
                    "reference_gene_id": gene,
                    "positive_unit_count": len(positive_units),
                    "not_called_unit_count": len(unknown_units),
                    "not_called_units": ";".join(
                        unit for unit in unit_order if unit in unknown_units
                    ),
                }
            )
        elif not positive_units:
            pattern = "no_loss"
            events = []
        else:
            events = maximal_loss_nodes(root, positive_units)
            if not events:
                raise ScaffoldError(f"{gene}: positive units produced no scaffold event")
            if len(events) > 1:
                pattern = "repeated_independent_events"
            elif events[0].node_type == "unit_terminal":
                pattern = "single_terminal_event"
            else:
                pattern = "single_internal_event"
        pattern_counts[pattern] += 1
        pattern_rows.append(
            {
                "reference_gene_id": gene,
                "scaffold_pattern": pattern,
                "positive_unit_count": len(positive_units),
                "event_count": len(events),
                "event_node_ids": ";".join(node.node_id for node in events),
            }
        )
        for node in events:
            descendants = node_descendants[node.node_id]
            supported_units = strict_gene_units.get(gene, set())
            event_rows.append(
                {
                    "reference_gene_id": gene,
                    "node_id": node.node_id,
                    "node_type": node.node_type,
                    "descendant_unit_count": len(descendants),
                    "descendant_units": ";".join(descendants),
                    "confirmed_decayed_supported_any_descendant": any(
                        unit in supported_units for unit in descendants
                    ),
                    "confirmed_decayed_supported_all_descendants": all(
                        unit in supported_units for unit in descendants
                    ),
                }
            )
    if sum(pattern_counts.values()) != args.expected_reference_genes:
        raise ScaffoldError("scaffold patterns do not close to reference genes")

    leaf_order = [leaf.name for leaf in root.leaves()]
    leaf_plot_order = {unit: order for order, unit in enumerate(leaf_order)}
    node_rows: list[dict[str, object]] = []
    for node in walk(root):
        descendants = node_descendants[node.node_id]
        node_rows.append(
            {
                "node_id": node.node_id,
                "parent_node_id": node.parent.node_id if node.parent else "",
                "node_type": node.node_type,
                "node_name": node.name,
                "depth": node_depth(node),
                "descendant_unit_count": len(descendants),
                "descendant_units": ";".join(descendants),
                "minimum_leaf_plot_order": min(
                    leaf_plot_order[unit] for unit in descendants
                ),
                "maximum_leaf_plot_order": max(
                    leaf_plot_order[unit] for unit in descendants
                ),
            }
        )
    node_frame = pd.DataFrame(node_rows)
    event_frame = pd.DataFrame(event_rows)
    pattern_frame = pd.DataFrame(pattern_rows)
    ambiguous_frame = pd.DataFrame(ambiguous_rows)
    branch_counts = Counter(event_frame["node_id"]) if not event_frame.empty else Counter()
    branch_strict_any = (
        event_frame.loc[
            event_frame["confirmed_decayed_supported_any_descendant"] == True,  # noqa: E712
            "node_id",
        ].value_counts().to_dict()
        if not event_frame.empty
        else {}
    )
    branch_strict_all = (
        event_frame.loc[
            event_frame["confirmed_decayed_supported_all_descendants"] == True,  # noqa: E712
            "node_id",
        ].value_counts().to_dict()
        if not event_frame.empty
        else {}
    )
    branch_frame = node_frame.copy()
    branch_frame["loss_event_gene_count"] = branch_frame["node_id"].map(
        lambda value: branch_counts[value]
    )
    branch_frame["confirmed_decayed_supported_any_descendant"] = branch_frame[
        "node_id"
    ].map(lambda value: int(branch_strict_any.get(value, 0)))
    branch_frame["confirmed_decayed_supported_all_descendants"] = branch_frame[
        "node_id"
    ].map(lambda value: int(branch_strict_all.get(value, 0)))
    pattern_summary = pd.DataFrame(
        [
            {
                "scaffold_pattern": key,
                "reference_gene_count": pattern_counts[key],
            }
            for key in (
                "no_loss",
                "single_terminal_event",
                "single_internal_event",
                "repeated_independent_events",
                "ambiguous_not_called",
            )
        ]
    )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        paths: list[Path] = []
        for filename, frame, compressed in (
            ("unit_loss_summary.tsv", unit_frame, False),
            ("unit_shared_nonshared_summary.tsv", shared_frame, False),
            ("unit_decayed_mechanism_summary.tsv", mechanism_frame, False),
            ("unit_decayed_cause_summary.tsv", cause_frame, False),
            ("loss_mechanism_crosswalk.tsv.gz", crosswalk, True),
            ("scaffold_nodes.tsv", node_frame, False),
            ("scaffold_branch_summary.tsv", branch_frame, False),
            ("scaffold_pattern_summary.tsv", pattern_summary, False),
            ("gene_scaffold_patterns.tsv.gz", pattern_frame, True),
            ("scaffold_loss_events.tsv.gz", event_frame, True),
            ("unplaced_not_called_genes.tsv", ambiguous_frame, False),
        ):
            path = staging / filename
            write_frame(path, frame, gzip_output=compressed)
            paths.append(path)
        tree_path = staging / "assembly_unit_scaffold.tre"
        tree_path.write_text(newick_record(root) + ";\n", encoding="utf-8")
        paths.append(tree_path)
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_UNIT_RESOLVED_ARTICLE_LOSS_SCAFFOLD",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "definitions": {
                "positive_loss": "manuscript_classification in {decayed, deleted}",
                "shared_loss": "positive in all assembly-unit terminals",
                "confirmed_decayed_types": (
                    "frameshift-only, in-frame-stop-only, and combined calls "
                    "intersected with manuscript_classification=decayed"
                ),
                "scaffold": (
                    "accepted 13-lineage topology expanded to parallel assembly-unit "
                    "terminals within multi-unit species; topology only, no dates or "
                    "branch-length interpretation"
                ),
                "event_placement": (
                    "maximal all-positive clades among genes with no not-called unit"
                ),
            },
            "counts": {
                "assembly_units": args.expected_units,
                "biological_lineages": args.expected_lineages,
                "reference_genes": args.expected_reference_genes,
                "matrix_rows": expected_rows,
                "positive_unit_gene_rows": args.expected_positive_rows,
                "shared_all_unit_genes": args.expected_shared_genes,
                "confirmed_decayed_unit_gene_rows": int(
                    mechanism_frame["confirmed_type_total"].sum()
                ),
                "scaffold_nodes": len(node_frame),
                "scaffold_terminal_nodes": int(
                    (node_frame["node_type"] == "unit_terminal").sum()
                ),
                "scaffold_event_rows": len(event_frame),
                "unplaced_not_called_genes": len(ambiguous_frame),
            },
            "confirmed_decayed_type_totals": expected_confirmed,
            "scaffold_pattern_counts": dict(pattern_counts),
            "inputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in inputs
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in paths
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
    except (
        ScaffoldError,
        OSError,
        csv.Error,
        UnicodeError,
        json.JSONDecodeError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
