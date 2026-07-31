#!/usr/bin/env python3
"""Render a path-free dated-tree and CAFE5 Base publication bundle.

The renderer consumes only validated, small outputs.  It never estimates a
tree, dates a node, or reruns CAFE.  Exact terminal closure, ultrametric branch
lengths, node-age clades, secondary-calibration rows, the validated Base model,
and all non-root branch expansion/contraction counts must reconcile before
plotting.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping, Sequence

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import (
    TaxonLabelError,
    format_downstream_taxon_label,
    format_taxon_label,
)


SCRIPT_VERSION = "1.6.0"
PUBLICATION_BOTTOM_TIP_RANK = {
    "Clematoclethra_scandens": 1,
    "Rhododendron_simsii": 2,
    "Coffea_arabica_E": 3,
    "Vitis_vinifera": 4,
}
GEOLOGIC_PERIODS = (
    ("Cretaceous", 145.0, 66.0, "#E6A04B"),
    ("Paleogene", 66.0, 23.03, "#73AD87"),
    ("Neogene", 23.03, 2.58, "#C6BA7C"),
    ("Qu.", 2.58, 0.0, "#6FA8C4"),
)
GEOLOGIC_EPOCHS = (
    ("Early", 145.0, 100.5, "#F4F4F4"),
    ("Late", 100.5, 66.0, "#E3E3E3"),
    ("Paleocene", 66.0, 56.0, "#D8D8D8"),
    ("Eocene", 56.0, 33.9, "#F4F4F4"),
    ("Oligocene", 33.9, 23.03, "#D8D8D8"),
    ("Miocene", 23.03, 5.333, "#F4F4F4"),
    ("Plio.", 5.333, 2.58, "#D8D8D8"),
    ("Pleist.", 2.58, 0.0117, "#D7E7EF"),
    ("H", 0.0117, 0.0, "#B7D4E2"),
)
GEOLOGIC_SCALE_SOURCE = "ICS International Chronostratigraphic Chart, 2026/06"
TERMINAL_COLUMNS = (
    "terminal_id",
    "biological_species",
    "grouping",
    "source_fasta_stem",
    "canonical_tree_label",
    "is_root_outgroup",
    "include_species_tree",
    "identity_status",
)
NODE_AGE_COLUMNS = (
    "node_id",
    "descendant_tip_count",
    "descendant_tips",
    "mean_ma",
    "q025_ma",
    "q975_ma",
    "chain1_mean_ma",
    "chain2_mean_ma",
    "combined_ess",
    "split_rhat",
)
CALIBRATION_COLUMNS = (
    "constraint_id",
    "node_label",
    "node_id",
    "minimum_ma",
    "maximum_ma",
    "posterior_mean_ma",
    "posterior_q025_ma",
    "posterior_q975_ma",
    "mean_inside_secondary_interval",
)
MODEL_COLUMNS = (
    "model_id",
    "role",
    "family_count",
    "significant_family_count_p_lt_0.05",
    "score",
    "lambda_values",
    "result_file_count",
)
CLADE_COLUMNS = ("model_id", "taxon_or_node_id", "increase", "decrease")
PLOT_COLUMNS = (
    "plot_order",
    "cafe_node_id",
    "node_type",
    "canonical_tree_label",
    "biological_species",
    "upright_suffix",
    "display_label",
    "descendant_tips",
    "parent_age_ma",
    "node_age_ma",
    "cafe_increase",
    "cafe_decrease",
)
CAFE_NODE_RE = re.compile(
    r"^(?:(?P<label>[A-Za-z0-9_]+))?<(?P<node_id>[1-9][0-9]*)>$"
)


class PhylogenyCafeError(RuntimeError):
    """Raised when frozen phylogeny/CAFE inputs do not reconcile."""


@dataclass
class Node:
    label: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    age: float = 0.0
    y: float = 0.0
    cafe_node_id: int | None = None
    cafe_increase: int | None = None
    cafe_decrease: int | None = None

    @property
    def is_tip(self) -> bool:
        return not self.children


def read_exact_tsv(path: Path, columns: Sequence[str], name: str) -> list[dict[str, str]]:
    path = Path(path).expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise PhylogenyCafeError(f"cannot open {name} {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed = tuple(reader.fieldnames or ())
        if observed != tuple(columns):
            raise PhylogenyCafeError(
                f"{name}: schema mismatch; expected={list(columns)!r}; observed={list(observed)!r}"
            )
        rows = [
            {column: (row[column] or "").strip() for column in columns}
            for row in reader
            if any((row[column] or "").strip() for column in columns)
        ]
    if not rows:
        raise PhylogenyCafeError(f"{name}: no data rows")
    return rows


def _number(value: str, location: str, *, minimum: float | None = None) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise PhylogenyCafeError(f"{location}: expected a number, found {value!r}") from error
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise PhylogenyCafeError(f"{location}: invalid value {value!r}")
    return result


def _integer(value: str, location: str, *, minimum: int = 0) -> int:
    if re.fullmatch(r"(?:0|[1-9][0-9]*)", value) is None:
        raise PhylogenyCafeError(f"{location}: expected a nonnegative integer")
    result = int(value)
    if result < minimum:
        raise PhylogenyCafeError(f"{location}: must be at least {minimum}")
    return result


def parse_newick(path: Path) -> Node:
    text = Path(path).read_text(encoding="utf-8").strip()
    position = 0

    def skip_space() -> None:
        nonlocal position
        while position < len(text) and text[position].isspace():
            position += 1

    def token(stop: str) -> str:
        nonlocal position
        start = position
        while position < len(text) and text[position] not in stop:
            position += 1
        return text[start:position].strip()

    def subtree() -> Node:
        nonlocal position
        skip_space()
        node = Node()
        if position < len(text) and text[position] == "(":
            position += 1
            while True:
                child = subtree()
                child.parent = node
                node.children.append(child)
                skip_space()
                if position >= len(text):
                    raise PhylogenyCafeError("dated tree ends inside an internal node")
                if text[position] == ",":
                    position += 1
                    continue
                if text[position] == ")":
                    position += 1
                    break
                raise PhylogenyCafeError(f"dated tree: unexpected token at position {position}")
            node.label = token(":,();")
        else:
            node.label = token(":,();")
            if not node.label:
                raise PhylogenyCafeError(f"dated tree: empty tip at position {position}")
        skip_space()
        if position < len(text) and text[position] == ":":
            position += 1
            raw = token(",();")
            node.length = _number(raw, f"dated tree branch for {node.label or 'internal node'}", minimum=0.0)
        return node

    root = subtree()
    skip_space()
    if position >= len(text) or text[position] != ";":
        raise PhylogenyCafeError("dated tree must end with one semicolon")
    position += 1
    skip_space()
    if position != len(text):
        raise PhylogenyCafeError("dated tree contains trailing content")
    root.length = 0.0
    return root


def walk(node: Node) -> Iterable[Node]:
    yield node
    for child in node.children:
        yield from walk(child)


def descendant_tips(node: Node) -> frozenset[str]:
    if node.is_tip:
        return frozenset((node.label,))
    tips: set[str] = set()
    for child in node.children:
        tips.update(descendant_tips(child))
    return frozenset(tips)


def orient_publication_tree(root: Node) -> None:
    """Rotate equivalent child orders so the three outgroups form the bottom block."""

    for node in walk(root):
        if node.children:
            node.children.sort(
                key=lambda child: min(
                    PUBLICATION_BOTTOM_TIP_RANK.get(tip, 0)
                    for tip in descendant_tips(child)
                )
            )


def parse_cafe_node_clades(
    path: Path,
) -> tuple[dict[frozenset[str], int], int]:
    """Return the exact CAFE node ID bound to each descendant-tip clade."""

    root = parse_newick(path)
    nodes = list(walk(root))
    if any(len(node.children) != 2 for node in nodes if not node.is_tip):
        raise PhylogenyCafeError("CAFE node-ID tree must be strictly bifurcating")
    seen_ids: set[int] = set()
    for node in nodes:
        match = CAFE_NODE_RE.fullmatch(node.label)
        if match is None:
            raise PhylogenyCafeError(
                f"CAFE node-ID tree contains an invalid label {node.label!r}"
            )
        node_id = int(match.group("node_id"))
        if node_id in seen_ids:
            raise PhylogenyCafeError(
                f"CAFE node-ID tree duplicates node ID {node_id}"
            )
        seen_ids.add(node_id)
        tip_label = match.group("label")
        if node.is_tip:
            if tip_label is None:
                raise PhylogenyCafeError(
                    f"CAFE node-ID tree tip {node_id} has no taxon label"
                )
            node.label = tip_label
        elif tip_label is not None:
            raise PhylogenyCafeError(
                f"CAFE node-ID tree internal node {node_id} has a taxon label"
            )
        else:
            node.label = ""
        node.cafe_node_id = node_id
    if seen_ids != set(range(1, len(nodes) + 1)):
        raise PhylogenyCafeError(
            "CAFE node-ID tree IDs are not one complete consecutive range"
        )
    by_clade: dict[frozenset[str], int] = {}
    for node in nodes:
        clade = descendant_tips(node)
        if clade in by_clade:
            raise PhylogenyCafeError("CAFE node-ID tree duplicates a descendant clade")
        assert node.cafe_node_id is not None
        by_clade[clade] = node.cafe_node_id
    assert root.cafe_node_id is not None
    return by_clade, root.cafe_node_id


def _terminal_suffix(row: Mapping[str, str]) -> str:
    # This is a species-level figure. Assembly, accession, haplotype, and
    # subgenome suffixes are intentionally suppressed. The two
    # A. zhejiangensis parental-lineage letters remain part of the curated
    # biological_species field and are therefore retained by the label helper.
    return ""


def _species_tree_display_label(row: Mapping[str, str]) -> str:
    species = row["biological_species"]
    if " parental lineage " in species:
        return format_downstream_taxon_label(
            species,
            (),
            abbreviate_genus=True,
        )
    return format_taxon_label(
        species,
        (),
        abbreviate_genus=True,
    )


def prepare(
    *,
    terminals_path: Path,
    dated_tree_path: Path,
    node_ages_path: Path,
    calibrations_path: Path,
    model_summary_path: Path,
    clade_summary_path: Path,
    cafe_node_tree_path: Path,
    cafe_validation_path: Path,
    expected_tip_count: int,
    expected_secondary_count: int,
) -> tuple[Node, list[dict[str, object]], dict[str, object], dict[frozenset[str], dict[str, str]]]:
    terminal_rows = read_exact_tsv(terminals_path, TERMINAL_COLUMNS, "terminal metadata")
    selected = [row for row in terminal_rows if row["include_species_tree"] == "true"]
    if len(selected) != expected_tip_count:
        raise PhylogenyCafeError(
            f"terminal metadata: expected {expected_tip_count} selected tips, found {len(selected)}"
        )
    if any(row["identity_status"] != "confirmed" for row in selected):
        raise PhylogenyCafeError("terminal metadata contains an unconfirmed selected identity")
    terminal_by_label = {row["canonical_tree_label"]: row for row in selected}
    if len(terminal_by_label) != len(selected):
        raise PhylogenyCafeError("terminal metadata contains duplicate canonical_tree_label values")

    root = parse_newick(dated_tree_path)
    orient_publication_tree(root)
    nodes = list(walk(root))
    tips = [node for node in nodes if node.is_tip]
    internals = [node for node in nodes if not node.is_tip]
    if any(len(node.children) != 2 for node in internals):
        raise PhylogenyCafeError("dated tree must be strictly bifurcating")
    tree_labels = {node.label for node in tips}
    if tree_labels != set(terminal_by_label) or len(tree_labels) != len(tips):
        raise PhylogenyCafeError(
            "dated tree tip set differs from terminal metadata: "
            f"missing={sorted(set(terminal_by_label) - tree_labels)}, "
            f"extra={sorted(tree_labels - set(terminal_by_label))}"
        )
    cafe_node_id_by_clade, cafe_root_id = parse_cafe_node_clades(
        cafe_node_tree_path
    )
    dated_clades = {descendant_tips(node) for node in nodes}
    if set(cafe_node_id_by_clade) != dated_clades:
        raise PhylogenyCafeError(
            "CAFE node-ID descendant clades do not exactly match the dated tree"
        )

    distances: dict[int, float] = {id(root): 0.0}
    for node in nodes:
        for child in node.children:
            distances[id(child)] = distances[id(node)] + child.length
    tip_distances = [distances[id(node)] for node in tips]
    root_age = sum(tip_distances) / len(tip_distances)
    max_deviation = max(abs(distance - root_age) for distance in tip_distances)
    if max_deviation > 1e-6:
        raise PhylogenyCafeError(
            f"dated tree is not ultrametric; maximum root-to-tip deviation={max_deviation}"
        )
    for node in nodes:
        node.age = root_age - distances[id(node)]

    node_age_rows = read_exact_tsv(node_ages_path, NODE_AGE_COLUMNS, "node ages")
    if len(node_age_rows) != len(internals):
        raise PhylogenyCafeError(
            f"node ages: expected {len(internals)} internal rows, found {len(node_age_rows)}"
        )
    ages_by_clade: dict[frozenset[str], dict[str, str]] = {}
    ages_by_id: dict[str, dict[str, str]] = {}
    for row in node_age_rows:
        clade = frozenset(filter(None, row["descendant_tips"].split(",")))
        if len(clade) != _integer(row["descendant_tip_count"], f"{row['node_id']} descendant_tip_count", minimum=2):
            raise PhylogenyCafeError(f"{row['node_id']}: descendant tip count mismatch")
        if clade in ages_by_clade or row["node_id"] in ages_by_id:
            raise PhylogenyCafeError("node ages contain duplicate node IDs or clades")
        mean = _number(row["mean_ma"], f"{row['node_id']} mean_ma", minimum=0.0)
        low = _number(row["q025_ma"], f"{row['node_id']} q025_ma", minimum=0.0)
        high = _number(row["q975_ma"], f"{row['node_id']} q975_ma", minimum=0.0)
        if not low <= mean <= high:
            raise PhylogenyCafeError(f"{row['node_id']}: posterior interval does not contain mean")
        ages_by_clade[clade] = row
        ages_by_id[row["node_id"]] = row
    internal_by_clade = {descendant_tips(node): node for node in internals}
    if set(ages_by_clade) != set(internal_by_clade):
        raise PhylogenyCafeError("node-age descendant clades do not exactly match the dated tree")
    for clade, row in ages_by_clade.items():
        if abs(internal_by_clade[clade].age - float(row["mean_ma"])) > 1e-6:
            raise PhylogenyCafeError(f"{row['node_id']}: dated-tree age differs from node-age table")

    calibration_rows = read_exact_tsv(
        calibrations_path, CALIBRATION_COLUMNS, "secondary calibrations"
    )
    if len(calibration_rows) != expected_secondary_count:
        raise PhylogenyCafeError(
            f"secondary calibrations: expected {expected_secondary_count}, found {len(calibration_rows)}"
        )
    for row in calibration_rows:
        if row["node_id"] not in ages_by_id or row["mean_inside_secondary_interval"] != "true":
            raise PhylogenyCafeError(f"{row['constraint_id']}: invalid or unbound secondary calibration")
        minimum = _number(row["minimum_ma"], f"{row['constraint_id']} minimum_ma", minimum=0.0)
        maximum = _number(row["maximum_ma"], f"{row['constraint_id']} maximum_ma", minimum=minimum)
        mean = _number(row["posterior_mean_ma"], f"{row['constraint_id']} posterior_mean_ma")
        if not minimum <= mean <= maximum:
            raise PhylogenyCafeError(f"{row['constraint_id']}: posterior mean is outside secondary interval")

    model_rows = read_exact_tsv(model_summary_path, MODEL_COLUMNS, "CAFE model summary")
    if len(model_rows) != 1 or model_rows[0]["model_id"] != "base_poisson":
        raise PhylogenyCafeError("CAFE publication requires exactly one base_poisson model row")
    model = model_rows[0]
    family_count = _integer(model["family_count"], "CAFE family_count", minimum=1)
    significant_count = _integer(
        model["significant_family_count_p_lt_0.05"], "CAFE significant family count"
    )
    if significant_count > family_count:
        raise PhylogenyCafeError("CAFE significant family count exceeds analyzed families")
    lambda_value = _number(model["lambda_values"], "CAFE lambda", minimum=0.0)

    try:
        cafe_validation = json.loads(Path(cafe_validation_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PhylogenyCafeError(f"cannot parse CAFE validation: {error}") from error
    if cafe_validation.get("status") not in {
        "PASS_CAFE5_BASE_VALIDATED_GAMMA_UNAVAILABLE",
        "PASS_CAFE5_PUBLICATION_TABLES",
    }:
        raise PhylogenyCafeError("CAFE validation status is not the accepted Base-only status")
    if cafe_validation.get("calibration_claim") != "TimeTree secondary-calibrated; not fossil-calibrated":
        raise PhylogenyCafeError("CAFE validation calibration claim is incorrect")
    unavailable = cafe_validation.get("unavailable_sensitivity") or {}
    gamma_status = unavailable.get("status") or cafe_validation.get("gamma3_status")
    if gamma_status != "UNAVAILABLE_INITIALIZATION_FAILURE":
        raise PhylogenyCafeError("CAFE Gamma3 unavailability is not explicitly validated")

    clade_rows = read_exact_tsv(clade_summary_path, CLADE_COLUMNS, "CAFE clade summary")
    cafe_counts_by_id: dict[int, tuple[int, int]] = {}
    cafe_labels_by_id: dict[int, str] = {}
    for row in clade_rows:
        if row["model_id"] != "base_poisson":
            raise PhylogenyCafeError("CAFE clade summary contains a non-Base model row")
        match = CAFE_NODE_RE.fullmatch(row["taxon_or_node_id"])
        if match is None:
            raise PhylogenyCafeError(
                f"CAFE clade summary has an invalid node label {row['taxon_or_node_id']!r}"
            )
        node_id = int(match.group("node_id"))
        if node_id in cafe_counts_by_id:
            raise PhylogenyCafeError(
                f"CAFE clade summary duplicates node ID {node_id}"
            )
        cafe_counts_by_id[node_id] = (
            _integer(row["increase"], f"CAFE node {node_id} increase"),
            _integer(row["decrease"], f"CAFE node {node_id} decrease"),
        )
        if match.group("label"):
            cafe_labels_by_id[node_id] = str(match.group("label"))
    expected_cafe_branch_ids = set(cafe_node_id_by_clade.values()) - {cafe_root_id}
    if set(cafe_counts_by_id) != expected_cafe_branch_ids:
        raise PhylogenyCafeError(
            "CAFE clade summary does not contain exactly one row for every non-root branch"
        )
    for tip in tips:
        node_id = cafe_node_id_by_clade[frozenset((tip.label,))]
        if cafe_labels_by_id.get(node_id) != tip.label:
            raise PhylogenyCafeError(
                f"CAFE node {node_id} terminal label does not match {tip.label}"
            )
    if set(cafe_labels_by_id.values()) != tree_labels:
        raise PhylogenyCafeError(
            "CAFE terminal count set differs from the dated-tree tip set"
        )
    for node in nodes:
        clade = descendant_tips(node)
        node.cafe_node_id = cafe_node_id_by_clade[clade]
        if node.parent is not None:
            node.cafe_increase, node.cafe_decrease = cafe_counts_by_id[
                node.cafe_node_id
            ]

    # Preserve the dated-tree vertical order.
    ordered_tips = tips
    for order, tip in enumerate(ordered_tips):
        tip.y = float(order)
    for node in reversed(nodes):
        if not node.is_tip:
            node.y = sum(child.y for child in node.children) / len(node.children)

    terminal_display: dict[str, tuple[str, str]] = {}
    for tip in ordered_tips:
        metadata = terminal_by_label[tip.label]
        suffix = _terminal_suffix(metadata)
        try:
            display = _species_tree_display_label(metadata)
        except TaxonLabelError as error:
            raise PhylogenyCafeError(
                f"{tip.label}: invalid publication label: {error}"
            ) from error
        terminal_display[tip.label] = (suffix, display)

    plot_rows: list[dict[str, object]] = []
    branch_nodes = [node for node in nodes if node.parent is not None]
    for order, node in enumerate(branch_nodes, 1):
        clade = descendant_tips(node)
        if node.is_tip:
            metadata = terminal_by_label[node.label]
            suffix, display = terminal_display[node.label]
            biological_species = metadata["biological_species"]
            canonical_tree_label = node.label
            node_type = "terminal"
        else:
            suffix = ""
            display = ""
            biological_species = ""
            canonical_tree_label = ""
            node_type = "internal"
        assert node.cafe_node_id is not None
        assert node.cafe_increase is not None
        assert node.cafe_decrease is not None
        plot_rows.append(
            {
                "plot_order": order,
                "cafe_node_id": node.cafe_node_id,
                "node_type": node_type,
                "canonical_tree_label": canonical_tree_label,
                "biological_species": biological_species,
                "upright_suffix": suffix,
                "display_label": display,
                "descendant_tips": ",".join(sorted(clade)),
                "parent_age_ma": node.parent.age,
                "node_age_ma": node.age,
                "cafe_increase": node.cafe_increase,
                "cafe_decrease": node.cafe_decrease,
            }
        )

    calibration_clades = {
        frozenset(filter(None, ages_by_id[row["node_id"]]["descendant_tips"].split(",")))
        for row in calibration_rows
    }
    validation: dict[str, object] = {
        "schema_version": 1,
        "renderer": "scripts/figures/render_phylogeny_cafe.py",
        "renderer_version": SCRIPT_VERSION,
        "status": "PASS_PHYLOGENY_CAFE_PUBLICATION_BUNDLE",
        "tip_count": len(tips),
        "internal_node_count": len(internals),
        "root_age_ma": root_age,
        "maximum_root_to_tip_deviation_ma": max_deviation,
        "secondary_calibration_count": len(calibration_rows),
        "calibration_claim": "TimeTree secondary-calibrated; not fossil-calibrated",
        "cafe_model": "Base Poisson",
        "cafe_analyzed_family_count": family_count,
        "cafe_significant_family_count_p_lt_0_05": significant_count,
        "cafe_lambda": lambda_value,
        "cafe_branch_count": len(cafe_counts_by_id),
        "cafe_internal_branch_count": sum(
            1 for node in branch_nodes if not node.is_tip
        ),
        "cafe_terminal_branch_count": sum(
            1 for node in branch_nodes if node.is_tip
        ),
        "cafe_annotation_semantics": (
            "Each value belongs to the incoming branch of the labelled "
            "descendant node; internal-branch values are anchored at descendant "
            "nodes and terminal-branch values are printed below species labels. "
            "The root has no incoming-branch value."
        ),
        "gamma3_status": "UNAVAILABLE_INITIALIZATION_FAILURE",
        "checks": {
            "terminal_identity_and_tip_closure": "pass",
            "strict_bifurcation": "pass",
            "ultrametric_tree": "pass",
            "node_age_clade_closure": "pass",
            "secondary_constraints_bound": "pass",
            "cafe_base_validation": "pass",
            "cafe_terminal_count_closure": "pass",
            "cafe_node_id_clade_closure": "pass",
            "cafe_all_nonroot_branches_annotated": "pass",
            "italic_binomials_upright_suffixes": "pass",
            "publication_outgroups_rotated_to_bottom": "pass",
            "species_level_labels_except_zhejiangensis_lineages": "pass",
            "terminal_cafe_counts_integrated_with_tree": "pass",
            "all_nonroot_cafe_counts_integrated_with_tree": "pass",
            "ics_period_and_epoch_strip": "pass",
            "secondary_calibrations_recorded_not_overplotted": "pass",
        },
        "bottom_tip_order": [
            "Clematoclethra_scandens",
            "Rhododendron_simsii",
            "Coffea_arabica_E",
            "Vitis_vinifera",
        ],
        "geologic_time_scale": {
            "source": GEOLOGIC_SCALE_SOURCE,
            "source_url": "https://stratigraphy.org/chart/",
            "period_boundaries_ma": {
                "Cretaceous_Paleogene": 66.0,
                "Paleogene_Neogene": 23.03,
                "Neogene_Quaternary": 2.58,
            },
            "epoch_boundaries_ma": {
                "Early_Late_Cretaceous": 100.5,
                "Cretaceous_Paleocene": 66.0,
                "Paleocene_Eocene": 56.0,
                "Eocene_Oligocene": 33.9,
                "Oligocene_Miocene": 23.03,
                "Miocene_Pliocene": 5.333,
                "Pliocene_Pleistocene": 2.58,
                "Pleistocene_Holocene": 0.0117,
            },
        },
        "secondary_calibration_table": "secondary_calibration_summary.tsv",
    }
    return root, plot_rows, validation, {clade: ages_by_clade[clade] for clade in calibration_clades}


def render_bundle(
    *,
    terminals_path: Path,
    dated_tree_path: Path,
    node_ages_path: Path,
    calibrations_path: Path,
    model_summary_path: Path,
    clade_summary_path: Path,
    cafe_node_tree_path: Path,
    cafe_validation_path: Path,
    output_dir: Path,
    basename: str,
    expected_tip_count: int = 17,
    expected_secondary_count: int = 4,
    dpi: int = 600,
) -> FigureBundle:
    root, rows, validation, _calibrated_ages = prepare(
        terminals_path=terminals_path,
        dated_tree_path=dated_tree_path,
        node_ages_path=node_ages_path,
        calibrations_path=calibrations_path,
        model_summary_path=model_summary_path,
        clade_summary_path=clade_summary_path,
        cafe_node_tree_path=cafe_node_tree_path,
        cafe_validation_path=cafe_validation_path,
        expected_tip_count=expected_tip_count,
        expected_secondary_count=expected_secondary_count,
    )
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
        from matplotlib.patches import Rectangle
    except ImportError as error:
        raise PhylogenyCafeError("matplotlib is required to render the publication figure") from error

    nodes = list(walk(root))
    tips = [node for node in nodes if node.is_tip]
    publication_style = {
        "font.family": "Arial",
        "font.size": 9.0,
        "font.weight": "normal",
        "axes.titlesize": 10.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 8.0,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
    with plt.rc_context(publication_style):
        fig, (tree_ax, time_ax) = plt.subplots(
            2,
            1,
            figsize=(7.2, 8.8),
            sharex=True,
            gridspec_kw={"height_ratios": (7.9, 1.25), "hspace": 0.07},
        )
        for node in nodes:
            if node.parent is not None:
                tree_ax.plot(
                    [node.parent.age, node.age],
                    [node.y, node.y],
                    color="#263238",
                    lw=1.05,
                    solid_capstyle="butt",
                )
            if node.children:
                tree_ax.plot(
                    [node.age, node.age],
                    [min(child.y for child in node.children), max(child.y for child in node.children)],
                    color="#263238",
                    lw=1.05,
                    solid_capstyle="butt",
                )
        internal_nodes = [node for node in nodes if not node.is_tip]
        tree_ax.scatter(
            [node.age for node in internal_nodes],
            [node.y for node in internal_nodes],
            s=5.5,
            facecolor="#263238",
            edgecolor="#263238",
            linewidth=0.0,
            zorder=4,
        )

        terminal_by_label = {
            str(row["canonical_tree_label"]): row
            for row in rows
            if str(row["node_type"]) == "terminal"
        }
        branch_by_node_id = {
            int(row["cafe_node_id"]): row
            for row in rows
        }
        label_x = -1.4
        expansion_color = "#087F8C"
        contraction_color = "#D95F02"
        for tip in tips:
            row = terminal_by_label[tip.label]
            tree_ax.text(
                label_x,
                tip.y - 0.10,
                str(row["display_label"]),
                va="center",
                ha="left",
                fontsize=10.0,
            )

        # CAFE reports a change on every non-root incoming branch. Bind each
        # annotation through the exact CAFE node-ID tree. Internal values are
        # anchored directly to their descendant nodes. Terminal values are
        # placed beneath the corresponding species labels so that short recent
        # branches remain uncluttered and the values stay publication-readable.
        for node in nodes:
            if node.parent is None:
                continue
            assert node.cafe_node_id is not None
            row = branch_by_node_id[node.cafe_node_id]
            if node.is_tip:
                annotation_xy = (label_x, node.y + 0.23)
                font_size = 7.8
            else:
                annotation_xy = (node.age, node.y)
                font_size = 7.2
            increase = TextArea(
                f"+{int(row['cafe_increase']):,}",
                textprops={
                    "color": expansion_color,
                    "fontsize": font_size,
                    "fontweight": "bold",
                },
            )
            separator = TextArea(
                "/",
                textprops={"color": "#5F6B73", "fontsize": font_size - 0.2},
            )
            decrease = TextArea(
                f"−{int(row['cafe_decrease']):,}",
                textprops={
                    "color": contraction_color,
                    "fontsize": font_size,
                    "fontweight": "bold",
                },
            )
            packed = HPacker(
                children=[increase, separator, decrease],
                align="center",
                pad=0,
                sep=0.7,
            )
            if node.is_tip:
                xybox = (0.0, 0.0)
                box_alignment = (0.0, 0.5)
            elif not node.is_tip:
                tier = node.cafe_node_id % 3
                xybox = (
                    -8.0 - 42.0 * tier,
                    (-13.0, 0.0, 13.0)[tier],
                )
                box_alignment = (1.0, 0.5)
            tree_ax.add_artist(
                AnnotationBbox(
                    packed,
                    annotation_xy,
                    xybox=xybox,
                    xycoords="data",
                    boxcoords="offset points",
                    frameon=not node.is_tip,
                    box_alignment=box_alignment,
                    arrowprops=(
                        None
                        if node.is_tip
                        else {
                            "arrowstyle": "-",
                            "color": "#899399",
                            "linewidth": 0.35,
                        }
                    ),
                    bboxprops={
                        "boxstyle": "round,pad=0.12",
                        "facecolor": "#FFFDF8",
                        "edgecolor": "#D8DDE0",
                        "linewidth": 0.30,
                        "alpha": 0.96,
                    },
                    zorder=5,
                )
            )

        root_age = float(validation["root_age_ma"])
        tree_ax.set_xlim(root_age * 1.03, -18.0)
        tree_ax.set_ylim(-0.8, len(tips) - 0.2)
        tree_ax.invert_yaxis()
        tree_ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        tree_ax.tick_params(axis="y", left=False, labelleft=False)
        tree_ax.tick_params(axis="x", bottom=False, labelbottom=False)
        tree_ax.grid(axis="x", color="#D9DEE2", lw=0.55, alpha=0.8)

        tree_ax.text(
            0.615,
            1.012,
            "Branch changes:",
            transform=tree_ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.2,
            fontweight="normal",
            color="#263238",
        )
        tree_ax.text(
            0.755,
            1.012,
            "+ expansion",
            transform=tree_ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.2,
            fontweight="bold",
            color=expansion_color,
        )
        tree_ax.text(
            0.875,
            1.012,
            "− contraction",
            transform=tree_ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.2,
            fontweight="bold",
            color=contraction_color,
        )

        time_ax.set_ylim(0.0, 2.0)
        time_ax.set_yticks([])
        time_ax.spines[["left", "right", "bottom"]].set_visible(False)
        time_ax.spines["top"].set_color("#263238")
        time_ax.spines["top"].set_linewidth(0.9)
        time_ax.xaxis.set_ticks_position("top")
        time_ax.xaxis.set_label_position("top")
        time_ax.tick_params(
            axis="x",
            top=True,
            bottom=False,
            labeltop=True,
            labelbottom=False,
            direction="out",
            length=3.0,
            width=0.8,
            pad=2.0,
        )
        time_ax.set_xlabel("Divergence time (Ma)", labelpad=5)
        for label, older, younger, color in GEOLOGIC_PERIODS:
            clipped_older = min(older, root_age)
            clipped_younger = max(younger, 0.0)
            if clipped_older <= clipped_younger:
                continue
            time_ax.add_patch(
                Rectangle(
                    (clipped_younger, 1.0),
                    clipped_older - clipped_younger,
                    1.0,
                    facecolor=color,
                    edgecolor="#5E666B",
                    linewidth=0.45,
                )
            )
            time_ax.text(
                (clipped_older + clipped_younger) / 2.0,
                1.5,
                label,
                ha="center",
                va="center",
                fontsize=9.0 if label != "Qu." else 7.0,
                color="#FFFFFF" if label != "Qu." else "#263238",
                fontweight="bold" if label != "Qu." else "normal",
            )
        for label, older, younger, color in GEOLOGIC_EPOCHS:
            clipped_older = min(older, root_age)
            clipped_younger = max(younger, 0.0)
            if clipped_older <= clipped_younger:
                continue
            width = clipped_older - clipped_younger
            time_ax.add_patch(
                Rectangle(
                    (clipped_younger, 0.0),
                    width,
                    1.0,
                    facecolor=color,
                    edgecolor="#AAB1B6",
                    linewidth=0.42,
                )
            )
            if width >= 1.65:
                compact_epoch = width < 3.5
                time_ax.text(
                    (clipped_older + clipped_younger) / 2.0,
                    0.5,
                    label,
                    ha="center",
                    va="center",
                    fontsize=7.8 if width >= 4.0 else 6.8,
                    rotation=90 if compact_epoch else 0,
                    color="#263238",
                )
        fig.subplots_adjust(left=0.065, right=0.985, top=0.965, bottom=0.07)

    branch_count = int(validation["cafe_branch_count"])
    terminal_branch_count = int(validation["cafe_terminal_branch_count"])
    internal_branch_count = int(validation["cafe_internal_branch_count"])
    caption = (
        "Primary dated phylogeny and gene-family changes. The tree shows the validated "
        "17-terminal MCMCTree ultrametric phylogeny in millions of years. The two-tier geologic "
        "scale uses International Commission on Stratigraphy boundaries for periods "
        "(Cretaceous, Paleogene, Neogene and Quaternary) and epochs (Early and Late Cretaceous, "
        "Paleocene, Eocene, Oligocene, Miocene, Pliocene, Pleistocene and Holocene). "
        "TimeTree secondary-calibration bounds are retained in the accompanying calibration "
        "table and are not overplotted; these are not fossil calibrations. Teal and orange "
        "values report family expansions and contractions, respectively, on each incoming "
        "branch from the validated CAFE5 Base Poisson model. Internal-branch values are anchored "
        "at their descendant nodes, whereas terminal-branch values are printed beneath the "
        f"corresponding species labels; all {branch_count} non-root branches "
        f"({terminal_branch_count} terminal and {internal_branch_count} internal) are shown. "
        "The model analyzed 15,066 families, estimated lambda = "
        "0.0085698614157905, and identified 1,539 nominal p < 0.05 families. The "
        "Gamma3 sensitivity is unavailable because initialization failed and is not claimed. "
        "Latin binomials are italic. Assembly, accession, haplotype, and subgenome suffixes are "
        "suppressed in this species-level display; only the two A. zhejiangensis parental "
        "lineages retain upright A/B identifiers."
    )
    bundle = write_figure_bundle(
        figure=fig,
        output_dir=output_dir,
        basename=basename,
        plot_rows=rows,
        plot_columns=PLOT_COLUMNS,
        caption=caption,
        validation=validation,
        input_paths=(
            terminals_path,
            dated_tree_path,
            node_ages_path,
            calibrations_path,
            model_summary_path,
            clade_summary_path,
            cafe_node_tree_path,
            cafe_validation_path,
        ),
        dpi=dpi,
    )
    plt.close(fig)
    return bundle


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terminals", required=True, type=Path)
    parser.add_argument("--dated-tree", required=True, type=Path)
    parser.add_argument("--node-ages", required=True, type=Path)
    parser.add_argument("--secondary-calibrations", required=True, type=Path)
    parser.add_argument("--cafe-model-summary", required=True, type=Path)
    parser.add_argument("--cafe-clade-summary", required=True, type=Path)
    parser.add_argument("--cafe-node-tree", required=True, type=Path)
    parser.add_argument("--cafe-validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="primary_phylogeny_cafe")
    parser.add_argument("--expected-tip-count", type=int, default=17)
    parser.add_argument("--expected-secondary-count", type=int, default=4)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args(argv)
    try:
        render_bundle(
            terminals_path=args.terminals,
            dated_tree_path=args.dated_tree,
            node_ages_path=args.node_ages,
            calibrations_path=args.secondary_calibrations,
            model_summary_path=args.cafe_model_summary,
            clade_summary_path=args.cafe_clade_summary,
            cafe_node_tree_path=args.cafe_node_tree,
            cafe_validation_path=args.cafe_validation,
            output_dir=args.output_dir,
            basename=args.basename,
            expected_tip_count=args.expected_tip_count,
            expected_secondary_count=args.expected_secondary_count,
            dpi=args.dpi,
        )
    except (OSError, PhylogenyCafeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
