#!/usr/bin/env python3
"""Aggregate article-method unit calls into complete biological-species losses.

For a biological species represented by more than one assembly unit, a
reference gene is positive only when every constituent unit is classified as
``decayed`` or ``deleted`` under the historical manuscript rule.  Mixed
positive/retained profiles are retained as ``partial_homeolog_specific`` and
never enter the complete-loss foreground.  The complete species states are
then placed on the accepted 13-lineage Actinidia topology as maximal
all-positive clades.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import prepare_unit_loss_evidence_scaffold as tree_utils


ARTICLE_CLASSES = {"retained", "decayed", "deleted", "not_called_loss"}
POSITIVE = {"decayed", "deleted"}


class SpeciesScaffoldError(ValueError):
    """Raised when frozen inputs cannot support species-complete aggregation."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-matrix", required=True, type=Path)
    parser.add_argument("--resolved-background", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--tip-map", required=True, type=Path)
    parser.add_argument("--time-tree", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-lineages", type=int, default=13)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-background-genes", type=int, default=33998)
    parser.add_argument("--expected-matrix-sha256", default="")
    parser.add_argument("--expected-background-sha256", default="")
    parser.add_argument("--expected-time-tree-sha256", default="")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_background(path: Path, expected: int) -> list[str]:
    genes = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(genes) != expected or len(set(genes)) != expected:
        raise SpeciesScaffoldError(
            f"{path.name}: expected {expected} unique genes, found {len(set(genes))}"
        )
    return genes


def species_identifier(species: str) -> str:
    return tree_utils.slug(species)


def assign_species_node_ids(root: tree_utils.Node) -> dict[str, tree_utils.Node]:
    result: dict[str, tree_utils.Node] = {}
    for node in tree_utils.walk(root):
        descendants = sorted(leaf.name for leaf in node.leaves())
        if not node.children:
            node.node_type = "species_terminal"
            node.node_id = f"terminal__{node.name}"
        elif node is root:
            node.node_type = "backbone_internal"
            node.node_id = "internal__actinidia_all"
        else:
            node.node_type = "backbone_internal"
            token = hashlib.sha1(
                "\n".join(descendants).encode("utf-8")
            ).hexdigest()[:10]
            node.node_id = f"internal__{token}"
        if node.node_id in result:
            raise SpeciesScaffoldError(f"duplicate node ID: {node.node_id}")
        result[node.node_id] = node
    return result


def write_frame(path: Path, frame: pd.DataFrame, *, gz: bool = False) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "mtime": 0} if gz else None,
    )


def run(args: argparse.Namespace) -> None:
    inputs = {
        "manuscript_matrix": args.manuscript_matrix,
        "resolved_background": args.resolved_background,
        "unit_metadata": args.unit_metadata,
        "tip_map": args.tip_map,
        "time_tree": args.time_tree,
    }
    for role, path in inputs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise SpeciesScaffoldError(f"missing or empty {role}: {path}")
    if args.output_dir.exists():
        raise SpeciesScaffoldError(
            f"output directory already exists: {args.output_dir}"
        )
    hashes = {role: sha256(path) for role, path in inputs.items()}
    expected_hashes = {
        "manuscript_matrix": args.expected_matrix_sha256,
        "resolved_background": args.expected_background_sha256,
        "time_tree": args.expected_time_tree_sha256,
    }
    for role, expected in expected_hashes.items():
        if expected and hashes[role] != expected.lower():
            raise SpeciesScaffoldError(
                f"{role} SHA-256 mismatch: {hashes[role]} != {expected.lower()}"
            )

    metadata, species_units = tree_utils.read_metadata(
        args.unit_metadata,
        args.expected_units,
        args.expected_lineages,
    )
    unit_order = list(metadata["assembly_unit_id"])
    unit_set = set(unit_order)
    species_ids = {
        species: species_identifier(species) for species in species_units
    }
    if len(set(species_ids.values())) != args.expected_lineages:
        raise SpeciesScaffoldError("species identifiers are not unique")

    root_all = tree_utils.parse_newick(args.time_tree)
    root, tip_to_species = tree_utils.actinidia_subtree(
        root_all,
        args.tip_map,
        set(species_units),
    )
    for leaf in root.leaves():
        leaf.name = species_ids[tip_to_species[leaf.name]]
        leaf.node_type = "species_terminal"
    assign_species_node_ids(root)
    expected_species_set = set(species_ids.values())
    if (
        len(root.leaves()) != args.expected_lineages
        or {leaf.name for leaf in root.leaves()} != expected_species_set
    ):
        raise SpeciesScaffoldError("species tree does not close to 13 lineages")

    background_order = read_background(
        args.resolved_background, args.expected_background_genes
    )
    background = set(background_order)
    needed = [
        "reference_gene_id",
        "assembly_unit_id",
        "manuscript_classification",
        "manuscript_positive_loss",
    ]
    matrix = tree_utils.read_tsv(args.manuscript_matrix, usecols=needed)
    expected_rows = args.expected_units * args.expected_reference_genes
    if len(matrix) != expected_rows:
        raise SpeciesScaffoldError(
            f"matrix rows={len(matrix):,}; expected {expected_rows:,}"
        )
    if matrix.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
        raise SpeciesScaffoldError("matrix contains duplicate gene-unit rows")
    if set(matrix["assembly_unit_id"]) != unit_set:
        raise SpeciesScaffoldError("matrix unit set differs from metadata")
    if set(matrix["manuscript_classification"]) - ARTICLE_CLASSES:
        raise SpeciesScaffoldError("matrix contains an unsupported article class")
    declared = matrix["manuscript_positive_loss"].str.lower() == "true"
    derived = matrix["manuscript_classification"].isin(POSITIVE)
    if not declared.equals(derived):
        raise SpeciesScaffoldError("article class and positive flag disagree")

    matrix = matrix.loc[matrix["reference_gene_id"].isin(background)].copy()
    if len(matrix) != args.expected_units * args.expected_background_genes:
        raise SpeciesScaffoldError("resolved background does not form a full unit grid")
    if (matrix["manuscript_classification"] == "not_called_loss").any():
        raise SpeciesScaffoldError(
            "resolved background unexpectedly contains not-called unit states"
        )
    state_grid = matrix.pivot(
        index="reference_gene_id",
        columns="assembly_unit_id",
        values="manuscript_classification",
    )[unit_order]
    if set(state_grid.index) != background:
        raise SpeciesScaffoldError("matrix/background gene sets differ")
    state_grid = state_grid.loc[background_order]

    species_order = [leaf.name for leaf in root.leaves()]
    species_name_by_id = {
        identifier: species for species, identifier in species_ids.items()
    }
    species_state = pd.DataFrame(index=background_order)
    summary_rows: list[dict[str, object]] = []
    terminal_metadata_rows: list[dict[str, object]] = []
    unit_metadata_lookup = metadata.set_index("assembly_unit_id").to_dict("index")
    for species_id in species_order:
        species = species_name_by_id[species_id]
        units = species_units[species]
        subset = state_grid[units]
        positive_count = subset.isin(POSITIVE).sum(axis=1)
        retained_count = (subset == "retained").sum(axis=1)
        if not ((positive_count + retained_count) == len(units)).all():
            raise SpeciesScaffoldError(f"{species}: states do not close")
        states = pd.Series(
            "partial_homeolog_specific",
            index=state_grid.index,
            dtype="object",
        )
        states.loc[positive_count == len(units)] = "complete_loss"
        states.loc[positive_count == 0] = "retained"
        species_state[species_id] = states
        counts = Counter(states)
        if sum(counts.values()) != args.expected_background_genes:
            raise SpeciesScaffoldError(f"{species}: aggregate states do not close")
        summary_rows.append(
            {
                "species_id": species_id,
                "biological_species": species,
                "constituent_unit_count": len(units),
                "constituent_units": ";".join(units),
                "complete_loss": counts["complete_loss"],
                "partial_homeolog_specific": counts[
                    "partial_homeolog_specific"
                ],
                "retained": counts["retained"],
                "resolved_background": args.expected_background_genes,
                "complete_loss_rate": counts["complete_loss"]
                / args.expected_background_genes,
            }
        )
        terminal_metadata_rows.append(
            {
                "assembly_unit_id": species_id,
                "biological_species": species,
                "haplotype_or_subgenome": "",
                "assembly_scope": (
                    "biological_species_complete_loss_across_constituent_units"
                ),
                "include": "true",
                "constituent_unit_count": len(units),
                "constituent_units": ";".join(units),
            }
        )

    node_descendants = {
        node.node_id: tuple(leaf.name for leaf in node.leaves())
        for node in tree_utils.walk(root)
    }
    event_rows: list[dict[str, object]] = []
    pattern_rows: list[dict[str, object]] = []
    pattern_counts: Counter[str] = Counter()
    for gene, row in species_state.iterrows():
        positive_species = {
            species_id
            for species_id, state in row.items()
            if state == "complete_loss"
        }
        partial_species = {
            species_id
            for species_id, state in row.items()
            if state == "partial_homeolog_specific"
        }
        if not positive_species:
            pattern = (
                "partial_only_no_complete_loss"
                if partial_species
                else "no_complete_loss"
            )
            events: list[tree_utils.Node] = []
        else:
            events = tree_utils.maximal_loss_nodes(root, positive_species)
            if not events:
                raise SpeciesScaffoldError(f"{gene}: no event for positive species")
            if len(events) > 1:
                pattern = "repeated_independent_events"
            elif events[0].node_type == "species_terminal":
                pattern = "single_terminal_event"
            else:
                pattern = "single_internal_event"
        pattern_counts[pattern] += 1
        pattern_rows.append(
            {
                "reference_gene_id": gene,
                "species_complete_loss_count": len(positive_species),
                "species_partial_loss_count": len(partial_species),
                "event_pattern": pattern,
                "event_count": len(events),
                "event_node_ids": ";".join(node.node_id for node in events),
            }
        )
        for node in events:
            descendants = node_descendants[node.node_id]
            event_rows.append(
                {
                    "reference_gene_id": gene,
                    "node_id": node.node_id,
                    "node_type": node.node_type,
                    "descendant_species_count": len(descendants),
                    "descendant_units": ";".join(descendants),
                    "descendant_species_ids": ";".join(descendants),
                }
            )

    leaf_plot_order = {
        species_id: index for index, species_id in enumerate(species_order)
    }
    node_rows: list[dict[str, object]] = []
    for node in tree_utils.walk(root):
        descendants = node_descendants[node.node_id]
        node_rows.append(
            {
                "node_id": node.node_id,
                "parent_node_id": node.parent.node_id if node.parent else "",
                "node_type": node.node_type,
                "node_name": node.name,
                "depth": tree_utils.node_depth(node),
                "descendant_unit_count": len(descendants),
                "descendant_units": ";".join(descendants),
                "descendant_species_count": len(descendants),
                "descendant_species_ids": ";".join(descendants),
                "minimum_leaf_plot_order": min(
                    leaf_plot_order[item] for item in descendants
                ),
                "maximum_leaf_plot_order": max(
                    leaf_plot_order[item] for item in descendants
                ),
            }
        )
    node_frame = pd.DataFrame(node_rows)
    event_frame = pd.DataFrame(event_rows)
    branch_counts = (
        Counter(event_frame["node_id"]) if not event_frame.empty else Counter()
    )
    branch_frame = node_frame.copy()
    branch_frame["complete_loss_event_genes"] = branch_frame["node_id"].map(
        lambda value: branch_counts[value]
    )

    state_long = (
        species_state.rename_axis("reference_gene_id").reset_index()
        .melt(
            id_vars="reference_gene_id",
            var_name="species_id",
            value_name="species_loss_state",
        )
        .merge(
            pd.DataFrame(summary_rows)[
                [
                    "species_id",
                    "biological_species",
                    "constituent_unit_count",
                    "constituent_units",
                ]
            ],
            on="species_id",
            how="left",
            validate="many_to_one",
        )
    )
    summary_frame = pd.DataFrame(summary_rows)
    terminal_metadata = pd.DataFrame(terminal_metadata_rows)
    pattern_frame = pd.DataFrame(pattern_rows)
    pattern_summary = pd.DataFrame(
        [
            {
                "event_pattern": pattern,
                "reference_gene_count": count,
            }
            for pattern, count in sorted(pattern_counts.items())
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
        for filename, frame, gz in (
            ("species_complete_loss_summary.tsv", summary_frame, False),
            ("species_complete_loss_states.tsv.gz", state_long, True),
            ("species_terminal_metadata.tsv", terminal_metadata, False),
            ("species_scaffold_nodes.tsv", node_frame, False),
            ("species_scaffold_branch_summary.tsv", branch_frame, False),
            ("gene_species_scaffold_patterns.tsv.gz", pattern_frame, True),
            ("species_scaffold_pattern_summary.tsv", pattern_summary, False),
            ("species_complete_loss_events.tsv.gz", event_frame, True),
        ):
            path = staging / filename
            write_frame(path, frame, gz=gz)
            paths.append(path)
        background_path = staging / "resolved_background_gene_ids.txt"
        background_path.write_text(
            "\n".join(background_order) + "\n", encoding="utf-8"
        )
        paths.append(background_path)
        tree_path = staging / "species_complete_loss_scaffold.tre"
        tree_path.write_text(tree_utils.newick_record(root) + ";\n", encoding="utf-8")
        paths.append(tree_path)

        terminal_ids = set(terminal_metadata["assembly_unit_id"])
        terminal_event_memberships = int(
            event_frame.loc[
                event_frame["node_type"] == "species_terminal"
            ].shape[0]
        )
        if terminal_ids != {
            row["descendant_units"]
            for row in node_rows
            if row["node_type"] == "species_terminal"
        }:
            raise SpeciesScaffoldError("terminal metadata and tree tips differ")
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_SPECIES_COMPLETE_LOSS_SCAFFOLD",
            "created_at_utc": datetime.now(timezone.utc).isoformat(
                timespec="seconds"
            ),
            "definitions": {
                "unit_positive": (
                    "manuscript_classification in {decayed, deleted}"
                ),
                "species_complete_loss": (
                    "all constituent assembly units are unit-positive"
                ),
                "species_partial_homeolog_specific": (
                    "at least one constituent unit is unit-positive and at least "
                    "one is retained"
                ),
                "analysis_background": (
                    "reference genes with a resolved article-method state in all "
                    "23 assembly units"
                ),
                "event_placement": (
                    "maximal all-complete-loss clades on the accepted 13-lineage "
                    "biological-species topology"
                ),
            },
            "counts": {
                "assembly_units": args.expected_units,
                "biological_species": args.expected_lineages,
                "reference_genes_in_matrix": args.expected_reference_genes,
                "resolved_background_genes": args.expected_background_genes,
                "species_state_rows": len(state_long),
                "species_scaffold_nodes": len(node_frame),
                "complete_loss_event_rows": len(event_frame),
                "terminal_complete_loss_memberships": terminal_event_memberships,
            },
            "event_pattern_counts": dict(pattern_counts),
            "inputs": [
                {
                    "role": role,
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": hashes[role],
                }
                for role, path in inputs.items()
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
        staging.replace(args.output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (
        SpeciesScaffoldError,
        tree_utils.ScaffoldError,
        OSError,
        UnicodeError,
        pd.errors.ParserError,
    ) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
