#!/usr/bin/env python3
"""Prepare article-method unit trends and topology-aware loss foregrounds.

The article-comparable classification is deliberately simple: exact SynOrths
anchors are retained, candidates with a qualifying genome-wide tBLASTX hit are
decayed, and candidates without such a hit are deleted.  Decayed and deleted
are positive loss calls for this analysis.  The orthogonal refined-cause field
is never used to rewrite those historical classes.

The 23 assembly units remain the primary denominator for expression, copy, and
shared/non-shared descriptive trends.  Biological-species grouping is used
only for topology-aware complete-lineage loss: every assembly unit assigned to
a lineage must be positive before that lineage is treated as absent.  Mixed
positive/retained polyploid states are retained as partial lineage loss and are
not promoted to ancestral-branch events.
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


POSITIVE = {"decayed", "deleted"}
ARTICLE_CLASSES = {"retained", "decayed", "deleted", "not_called_loss"}
PLOIDY_LABELS = {"2x": "diploid", "4x": "tetraploid", "6x": "hexaploid"}
CLUSTER_MEMBER = re.compile(r">(.+?)\.\.\.")


class DownstreamError(ValueError):
    """Raised when frozen inputs cannot support the declared analysis."""


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None

    def leaves(self) -> list["Node"]:
        if not self.children:
            return [self]
        return [leaf for child in self.children for leaf in child.leaves()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-matrix", required=True, type=Path)
    parser.add_argument("--unit-ledger", required=True, type=Path)
    parser.add_argument("--species-map", required=True, type=Path)
    parser.add_argument("--tip-map", required=True, type=Path)
    parser.add_argument("--time-tree", required=True, type=Path)
    parser.add_argument("--reference-expression", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-lineages", type=int, default=13)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, compression="infer")


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise DownstreamError(f"{label} is missing columns: {', '.join(missing)}")


def read_units(path: Path, expected: int) -> pd.DataFrame:
    frame = read_tsv(path)
    require_columns(
        frame,
        {"sample_id", "species", "ploidy", "analysis_role", "input_scope"},
        "unit ledger",
    )
    frame = frame.loc[frame["analysis_role"] == "target_repertoire"].copy()
    if len(frame) != expected or frame["sample_id"].duplicated().any():
        raise DownstreamError(f"unit ledger must contain exactly {expected} unique targets")
    if set(frame["input_scope"]) != {"whole_genome"}:
        raise DownstreamError("every target unit must use whole_genome scope")
    unknown = sorted(set(frame["ploidy"]) - set(PLOIDY_LABELS))
    if unknown:
        raise DownstreamError(f"unsupported ploidy labels: {unknown}")
    frame["ploidy_class"] = frame["ploidy"].map(PLOIDY_LABELS)
    return frame[["sample_id", "species", "ploidy", "ploidy_class"]]


def read_species_map(path: Path, units: set[str], expected_lineages: int) -> dict[str, str]:
    frame = read_tsv(path)
    require_columns(frame, {"assembly_unit_id", "biological_species", "include"}, "species map")
    frame = frame.loc[frame["include"].str.lower() == "true"].copy()
    if frame["assembly_unit_id"].duplicated().any():
        raise DownstreamError("species map contains duplicate assembly units")
    mapping = dict(zip(frame["assembly_unit_id"], frame["biological_species"]))
    if set(mapping) != units:
        raise DownstreamError("species-map unit universe does not match the target ledger")
    if len(set(mapping.values())) != expected_lineages:
        raise DownstreamError(
            f"species map contains {len(set(mapping.values()))} lineages; expected {expected_lineages}"
        )
    return mapping


def tokenize_newick(text: str) -> list[str]:
    tokens = re.findall(r"\(|\)|,|:|;|[^\s(),:;]+", text.strip())
    if not tokens or tokens[-1] != ";":
        raise DownstreamError("time tree must be one semicolon-terminated Newick record")
    return tokens


def parse_newick(path: Path) -> Node:
    tokens = tokenize_newick(path.read_text(encoding="utf-8"))
    index = 0

    def subtree() -> Node:
        nonlocal index
        if index >= len(tokens):
            raise DownstreamError("unexpected end of Newick tree")
        if tokens[index] == "(":
            index += 1
            children = [subtree()]
            while index < len(tokens) and tokens[index] == ",":
                index += 1
                children.append(subtree())
            if index >= len(tokens) or tokens[index] != ")":
                raise DownstreamError("unbalanced Newick parentheses")
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
                raise DownstreamError(f"invalid Newick token {tokens[index]!r}")
            node = Node(name=tokens[index])
            index += 1
        if index < len(tokens) and tokens[index] == ":":
            index += 1
            if index >= len(tokens) or tokens[index] in {"(", ")", ",", ":", ";"}:
                raise DownstreamError("Newick branch length is missing")
            try:
                float(tokens[index])
            except ValueError as error:
                raise DownstreamError(f"invalid Newick branch length {tokens[index]!r}") from error
            index += 1
        return node

    root = subtree()
    if index != len(tokens) - 1:
        raise DownstreamError("unexpected tokens after Newick root")
    leaf_names = [leaf.name for leaf in root.leaves()]
    if any(not name for name in leaf_names) or len(leaf_names) != len(set(leaf_names)):
        raise DownstreamError("Newick tips must be nonempty and unique")
    return root


def read_tip_map(path: Path, species: set[str], root: Node) -> tuple[dict[str, str], Node]:
    frame = read_tsv(path)
    require_columns(frame, {"tree_tip", "biological_species", "include"}, "tip map")
    frame = frame.loc[frame["include"].str.lower() == "true"].copy()
    if frame["tree_tip"].duplicated().any() or frame["biological_species"].duplicated().any():
        raise DownstreamError("included tree tips and biological lineages must be one-to-one")
    tip_to_species = dict(zip(frame["tree_tip"], frame["biological_species"]))
    if set(tip_to_species.values()) != species:
        raise DownstreamError("tip-map biological lineages do not match the species map")
    tree_tips = {leaf.name for leaf in root.leaves()}
    if not set(tip_to_species).issubset(tree_tips):
        raise DownstreamError("tip map contains a terminal absent from the time tree")
    wanted = set(tip_to_species)

    def descendant_names(node: Node) -> set[str]:
        return {leaf.name for leaf in node.leaves()}

    candidates = [
        node for node in walk(root)
        if wanted.issubset(descendant_names(node))
    ]
    actinidia_root = min(candidates, key=lambda node: len(descendant_names(node)))
    if descendant_names(actinidia_root) != wanted:
        raise DownstreamError("included lineage tips do not form one exact Actinidia clade")
    return tip_to_species, actinidia_root


def walk(node: Node) -> list[Node]:
    return [node, *(descendant for child in node.children for descendant in walk(child))]


def branch_info(root: Node, tip_to_species: dict[str, str]) -> tuple[dict[int, str], dict[int, set[str]]]:
    identifiers: dict[int, str] = {}
    descendants: dict[int, set[str]] = {}
    for node in walk(root):
        taxa = {tip_to_species[leaf.name] for leaf in node.leaves()}
        descendants[id(node)] = taxa
        if len(taxa) == 1:
            species = next(iter(taxa))
            token = re.sub(r"[^a-z0-9]+", "_", species.lower()).strip("_")
            identifiers[id(node)] = f"terminal__{token}"
        elif node is root:
            identifiers[id(node)] = "internal__actinidia_all"
        else:
            digest = hashlib.sha1("\n".join(sorted(taxa)).encode("utf-8")).hexdigest()[:10]
            identifiers[id(node)] = f"internal__{digest}"
    return identifiers, descendants


def cluster_members(path: Path) -> set[str]:
    members: set[str] = set()
    current = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current = True
                continue
            match = CLUSTER_MEMBER.search(line)
            if not current or match is None:
                raise DownstreamError(f"malformed CD-HIT record at line {line_number}")
            gene = match.group(1)
            if gene in members:
                raise DownstreamError(f"duplicate CD-HIT gene: {gene}")
            members.add(gene)
    if not members:
        raise DownstreamError("no CD-HIT members parsed")
    return members


def classify_species(group: pd.Series) -> str:
    states = list(group)
    positive = sum(state in POSITIVE for state in states)
    retained = states.count("retained")
    not_called = states.count("not_called_loss")
    if positive == len(states):
        return "complete_loss"
    if retained and positive:
        return "partial_loss"
    if retained:
        return "retained"
    if not_called:
        return "unknown"
    raise DownstreamError(f"unhandled species state combination: {Counter(states)}")


def maximal_loss_nodes(root: Node, loss_species: set[str], descendants: dict[int, set[str]]) -> list[Node]:
    events: list[Node] = []

    def visit(node: Node, parent_all_loss: bool) -> None:
        taxa = descendants[id(node)]
        all_loss = bool(taxa) and taxa.issubset(loss_species)
        if all_loss and not parent_all_loss:
            events.append(node)
            return
        for child in node.children:
            visit(child, all_loss)

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
        args.manuscript_matrix, args.unit_ledger, args.species_map, args.tip_map,
        args.time_tree, args.reference_expression, args.clusters,
    ]
    missing = [str(path) for path in inputs if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise DownstreamError(f"missing or empty inputs: {missing}")
    if args.output_dir.exists():
        raise DownstreamError(f"output directory already exists: {args.output_dir}")

    units = read_units(args.unit_ledger, args.expected_units)
    unit_ids = set(units["sample_id"])
    species_map = read_species_map(args.species_map, unit_ids, args.expected_lineages)
    root_all = parse_newick(args.time_tree)
    tip_to_species, actinidia_root = read_tip_map(
        args.tip_map, set(species_map.values()), root_all
    )
    branch_ids, branch_descendants = branch_info(actinidia_root, tip_to_species)

    losses = read_tsv(args.manuscript_matrix)
    require_columns(
        losses,
        {
            "reference_gene_id", "assembly_unit_id", "manuscript_classification",
            "manuscript_positive_loss", "manuscript_rule", "refined_decayed_cause",
        },
        "manuscript matrix",
    )
    expected_rows = args.expected_units * args.expected_reference_genes
    if len(losses) != expected_rows:
        raise DownstreamError(f"manuscript matrix has {len(losses)} rows; expected {expected_rows}")
    if losses.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
        raise DownstreamError("manuscript matrix contains duplicate gene-unit rows")
    if set(losses["assembly_unit_id"]) != unit_ids:
        raise DownstreamError("manuscript matrix unit universe does not match the ledger")
    classes = set(losses["manuscript_classification"])
    if not classes.issubset(ARTICLE_CLASSES):
        raise DownstreamError(f"unexpected manuscript classes: {sorted(classes)}")
    expected_positive = losses["manuscript_classification"].isin(POSITIVE)
    observed_positive = losses["manuscript_positive_loss"].str.lower() == "true"
    if not expected_positive.equals(observed_positive):
        raise DownstreamError("manuscript positive flag is inconsistent with decayed+deleted")
    reference = sorted(losses["reference_gene_id"].unique())
    if len(reference) != args.expected_reference_genes:
        raise DownstreamError(
            f"manuscript matrix contains {len(reference)} genes; expected {args.expected_reference_genes}"
        )

    positive_counts = losses.groupby("reference_gene_id")["manuscript_positive_loss"].apply(
        lambda values: sum(value.lower() == "true" for value in values)
    )
    shared = set(positive_counts.loc[positive_counts == args.expected_units].index)
    if not shared:
        raise DownstreamError("no 23-unit shared article-method losses were found")
    nonshared = set(reference) - shared

    losses["species"] = losses["assembly_unit_id"].map(species_map)
    losses = losses.merge(
        units, left_on="assembly_unit_id", right_on="sample_id", validate="many_to_one",
        suffixes=("", "_ledger"),
    )
    if (losses["species"] != losses["species_ledger"]).any():
        raise DownstreamError("species labels disagree between unit ledger and aggregation map")

    class_counts = (
        losses.groupby(["assembly_unit_id", "manuscript_classification"]).size()
        .unstack(fill_value=0)
        .reindex(columns=["retained", "decayed", "deleted", "not_called_loss"], fill_value=0)
        .reset_index()
    )
    class_counts = class_counts.merge(
        units, left_on="assembly_unit_id", right_on="sample_id", validate="one_to_one"
    )
    class_counts["positive_loss"] = class_counts["decayed"] + class_counts["deleted"]
    class_counts["resolved_denominator"] = (
        class_counts["retained"] + class_counts["decayed"] + class_counts["deleted"]
    )
    class_counts["positive_loss_rate"] = (
        class_counts["positive_loss"] / class_counts["resolved_denominator"]
    )
    class_counts["shared_positive_loss"] = len(shared)
    class_counts["nonshared_positive_loss"] = class_counts["positive_loss"] - len(shared)
    class_counts["nonshared_positive_loss_rate"] = (
        class_counts["nonshared_positive_loss"]
        / (class_counts["resolved_denominator"] - len(shared))
    )
    class_counts = class_counts[
        [
            "assembly_unit_id", "species", "ploidy", "ploidy_class", "retained",
            "decayed", "deleted", "not_called_loss", "positive_loss",
            "resolved_denominator", "positive_loss_rate", "shared_positive_loss",
            "nonshared_positive_loss", "nonshared_positive_loss_rate",
        ]
    ].sort_values("assembly_unit_id")

    resolved_nonshared = losses.loc[
        losses["reference_gene_id"].isin(nonshared)
        & losses["manuscript_classification"].isin({"retained", "decayed", "deleted"})
    ].copy()
    resolved_nonshared = resolved_nonshared[
        [
            "reference_gene_id", "assembly_unit_id", "species", "ploidy_class", "ploidy",
            "manuscript_classification", "manuscript_rule", "refined_decayed_cause",
        ]
    ].rename(
        columns={
            "ploidy_class": "ploidy",
            "ploidy": "ploidy_level",
            "manuscript_classification": "classification",
            "manuscript_rule": "evidence_source",
            "refined_decayed_cause": "refined_cause",
        }
    ).sort_values(["assembly_unit_id", "reference_gene_id"])

    species_groups = defaultdict(list)
    for unit, species in species_map.items():
        species_groups[species].append(unit)
    species_rows: list[dict[str, object]] = []
    species_state_by_gene: dict[str, dict[str, str]] = defaultdict(dict)
    for (gene, species), group in losses.groupby(["reference_gene_id", "species"], sort=True):
        states = group["manuscript_classification"]
        status = classify_species(states)
        species_state_by_gene[gene][species] = status
        counts = Counter(states)
        species_rows.append(
            {
                "reference_gene_id": gene,
                "biological_species": species,
                "assembly_unit_count": len(group),
                "retained_unit_count": counts["retained"],
                "decayed_unit_count": counts["decayed"],
                "deleted_unit_count": counts["deleted"],
                "not_called_unit_count": counts["not_called_loss"],
                "species_gene_status": status,
            }
        )

    pattern_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    foregrounds: dict[str, set[str]] = defaultdict(set)
    foreground_meta: dict[str, dict[str, object]] = {}
    fully_resolved: set[str] = set()
    ordered_species = sorted(species_groups)
    for gene in reference:
        states = species_state_by_gene[gene]
        if set(states) != set(ordered_species):
            raise DownstreamError(f"incomplete species grid for {gene}")
        losses_complete = {species for species, state in states.items() if state == "complete_loss"}
        unknown = {species for species, state in states.items() if state == "unknown"}
        partial = {species for species, state in states.items() if state == "partial_loss"}
        exact = not unknown
        if exact:
            fully_resolved.add(gene)
        events = maximal_loss_nodes(actinidia_root, losses_complete, branch_descendants) if exact else []
        if unknown:
            pattern = "ambiguous_missing_data"
        elif not losses_complete and partial:
            pattern = "partial_lineage_loss_only"
        elif not losses_complete:
            pattern = "no_complete_lineage_loss"
        elif len(events) == 1 and len(branch_descendants[id(events[0])]) == 1:
            pattern = "single_terminal_branch_loss"
        elif len(events) == 1:
            pattern = "single_internal_branch_loss"
        else:
            pattern = "recurrent_independent_losses"
        pattern_rows.append(
            {
                "reference_gene_id": gene,
                "tree_pattern": pattern,
                "tree_placement_exact": str(exact).lower(),
                "complete_loss_lineage_count": len(losses_complete),
                "partial_loss_lineage_count": len(partial),
                "unknown_lineage_count": len(unknown),
                "minimum_loss_event_count": len(events) if exact else "",
                "complete_loss_lineages": ";".join(sorted(losses_complete)),
                "partial_loss_lineages": ";".join(sorted(partial)),
                "unknown_lineages": ";".join(sorted(unknown)),
            }
        )
        if partial and exact:
            foregrounds["category__partial_lineage_loss_any"].add(gene)
        if gene in shared:
            foregrounds["category__shared_all_23_units"].add(gene)
        if gene not in shared and positive_counts[gene] > 0:
            foregrounds["category__nonshared_any_unit_loss"].add(gene)
        if exact and events:
            foregrounds[f"category__{pattern}"].add(gene)
            for event_index, node in enumerate(events, 1):
                branch_id = branch_ids[id(node)]
                taxa = branch_descendants[id(node)]
                foreground_id = f"branch__{branch_id}"
                foregrounds[foreground_id].add(gene)
                foreground_meta[foreground_id] = {
                    "analysis_scope": "topology_inferred_complete_lineage_loss_branch",
                    "branch_id": branch_id,
                    "descendant_lineage_count": len(taxa),
                    "descendant_lineages": ";".join(sorted(taxa)),
                    "background_scope": "all_13_lineages_resolved_article_method",
                }
                event_rows.append(
                    {
                        "reference_gene_id": gene,
                        "event_index": event_index,
                        "branch_id": branch_id,
                        "branch_type": "terminal" if len(taxa) == 1 else "internal",
                        "descendant_lineage_count": len(taxa),
                        "descendant_lineages": ";".join(sorted(taxa)),
                        "gene_minimum_loss_event_count": len(events),
                    }
                )

    category_meta = {
        "category__shared_all_23_units": (
            "all_23_units_positive_article_method", "all_reference_genes"
        ),
        "category__nonshared_any_unit_loss": (
            "positive_in_at_least_one_but_not_all_23_units", "all_reference_genes"
        ),
        "category__partial_lineage_loss_any": (
            "mixed_positive_and_retained_units_within_at_least_one_lineage",
            "all_13_lineages_resolved_article_method",
        ),
        "category__single_terminal_branch_loss": (
            "one_exact_terminal_complete_lineage_loss_event",
            "all_13_lineages_resolved_article_method",
        ),
        "category__single_internal_branch_loss": (
            "one_exact_internal_complete_lineage_loss_event",
            "all_13_lineages_resolved_article_method",
        ),
        "category__recurrent_independent_losses": (
            "two_or_more_exact_complete_lineage_loss_events",
            "all_13_lineages_resolved_article_method",
        ),
    }
    for foreground_id, (definition, background_scope) in category_meta.items():
        foreground_meta[foreground_id] = {
            "analysis_scope": "topology_aware_article_method_category",
            "branch_id": "",
            "descendant_lineage_count": "",
            "descendant_lineages": "",
            "background_scope": background_scope,
            "definition": definition,
        }
    for foreground_id in foregrounds:
        foreground_meta.setdefault(
            foreground_id,
            {
                "analysis_scope": "topology_aware_article_method_category",
                "branch_id": "",
                "descendant_lineage_count": "",
                "descendant_lineages": "",
                "background_scope": "all_13_lineages_resolved_article_method",
                "definition": foreground_id.removeprefix("category__"),
            },
        )

    expression = read_tsv(args.reference_expression)
    require_columns(expression, {"reference_gene_id", "leaf_raw_count"}, "reference expression")
    if expression["reference_gene_id"].duplicated().any() or set(expression["reference_gene_id"]) != set(reference):
        raise DownstreamError("reference expression must contain the exact reference-gene universe")
    expression = expression.loc[expression["reference_gene_id"].isin(nonshared)].copy()
    expression["_order"] = expression["reference_gene_id"].map({gene: i for i, gene in enumerate(reference)})
    expression = expression.sort_values("_order").drop(columns="_order")
    members = cluster_members(args.clusters)
    if not members.issubset(set(reference)):
        raise DownstreamError("CD-HIT contains non-reference gene IDs")
    missing_copy = sorted(nonshared - members)

    parent = args.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.tmp.", dir=parent))
    try:
        shared_frame = pd.DataFrame(
            {
                "reference_gene_id": sorted(shared),
                "shared_positive_all_23_units": "true",
            }
        )
        species_frame = pd.DataFrame(species_rows).sort_values(
            ["reference_gene_id", "biological_species"]
        )
        patterns = pd.DataFrame(pattern_rows).sort_values("reference_gene_id")
        events = pd.DataFrame(event_rows)
        if events.empty:
            raise DownstreamError("no exact tree loss events were inferred")
        events = events.sort_values(["reference_gene_id", "event_index"])
        foreground_rows = [
            {"foreground_id": foreground_id, "reference_gene_id": gene}
            for foreground_id in sorted(foregrounds)
            for gene in sorted(foregrounds[foreground_id])
        ]
        metadata_rows = []
        for foreground_id in sorted(foregrounds):
            row = {"foreground_id": foreground_id, **foreground_meta[foreground_id]}
            row["foreground_gene_count"] = len(foregrounds[foreground_id])
            metadata_rows.append(row)

        write_frame(temporary / "unit_loss_summary.tsv", class_counts)
        write_frame(temporary / "shared_23_unit_genes.tsv", shared_frame)
        write_frame(
            temporary / "resolved_nonshared_unit_loss_table.tsv.gz",
            resolved_nonshared,
            gzip_output=True,
        )
        write_frame(temporary / "nonshared_reference_leaf_raw_counts.tsv", expression)
        (temporary / "nonshared_cdhit_missing_reference_ids.txt").write_text(
            "".join(f"{gene}\n" for gene in missing_copy), encoding="utf-8"
        )
        write_frame(
            temporary / "species_gene_states.tsv.gz", species_frame, gzip_output=True
        )
        write_frame(
            temporary / "gene_tree_patterns.tsv.gz", patterns, gzip_output=True
        )
        write_frame(
            temporary / "tree_loss_events.tsv.gz", events, gzip_output=True
        )
        write_frame(
            temporary / "tree_foreground_gene_ids.tsv.gz",
            pd.DataFrame(foreground_rows),
            gzip_output=True,
        )
        write_frame(
            temporary / "tree_foreground_metadata.tsv", pd.DataFrame(metadata_rows)
        )
        (temporary / "tree_background_gene_ids.txt").write_text(
            "".join(f"{gene}\n" for gene in sorted(fully_resolved)), encoding="utf-8"
        )
        output_names = sorted(path.name for path in temporary.iterdir() if path.is_file())
        checksum_rows = [
            {"file": name, "bytes": (temporary / name).stat().st_size, "sha256": sha256(temporary / name)}
            for name in output_names
        ]
        write_frame(temporary / "checksums.tsv", pd.DataFrame(checksum_rows))
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_MANUSCRIPT_METHOD_DOWNSTREAM_AND_TREE_PATTERNS",
            "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "definitions": {
                "unit_positive": "article-method decayed or deleted",
                "unit_denominator": "retained + decayed + deleted; not_called_loss excluded",
                "shared": "positive in all 23 assembly units",
                "complete_lineage_loss": "all assembly units assigned to a biological lineage are positive",
                "partial_lineage_loss": "at least one positive and at least one retained unit within a biological lineage",
                "tree_event": "minimum irreversible loss clades on the exact matching Actinidia topology; genes with unknown lineage states excluded",
                "refinement_policy": "refined causes never rewrite article-method retained/decayed/deleted classes",
            },
            "counts": {
                "reference_genes": len(reference),
                "assembly_units": len(unit_ids),
                "biological_lineages": len(ordered_species),
                "shared_23_unit_genes": len(shared),
                "nonshared_genes": len(nonshared),
                "genes_with_any_unit_positive": int((positive_counts > 0).sum()),
                "tree_fully_resolved_genes": len(fully_resolved),
                "tree_loss_event_rows": len(events),
                "tree_foregrounds": len(foregrounds),
                "cdhit_missing_nonshared_genes": len(missing_copy),
            },
            "tree_pattern_counts": dict(Counter(row["tree_pattern"] for row in pattern_rows)),
            "inputs": [
                {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                for path in inputs
            ],
            "outputs": checksum_rows,
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (DownstreamError, OSError, pd.errors.ParserError, UnicodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
