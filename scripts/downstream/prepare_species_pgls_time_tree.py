#!/usr/bin/env python3
"""Publish the exact biological-lineage PGLS subtree from a validated dated tree."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
from pathlib import Path

from Bio import Phylo


class TreePreparationError(RuntimeError):
    pass


CALIBRATION_CLAIM = "TimeTree secondary-calibrated; not fossil-calibrated"
PASS_CHECKS = {
    "dated_tree_manifest_pass",
    "rooted_by_accepted_topology",
    "exact_biological_species_tip_set",
    "strictly_bifurcating",
    "ultrametric",
    "finite_nonnegative_branch_lengths",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise TreePreparationError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise TreePreparationError(f"empty file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, str | int]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def validate_dated_directory(directory: Path) -> tuple[dict[str, object], Path]:
    if not directory.is_dir() or directory.is_symlink():
        raise TreePreparationError(f"invalid dating validation directory: {directory}")
    validation_path = regular(directory / "validation.json")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if (
        validation.get("status") != "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC"
        or validation.get("workflow")
        != "mcmctree_secondary_two_chain_validation_and_ultrametric_publication"
        or validation.get("calibration_claim") != CALIBRATION_CLAIM
    ):
        raise TreePreparationError("dated-tree validation is not the accepted TimeTree PASS")
    checksums_path = regular(directory / "checksums.tsv")
    with checksums_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise TreePreparationError("invalid dated-tree checksum header")
        rows = list(reader)
    checksums: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if not name or Path(name).name != name or name in checksums or len(digest) != 64:
            raise TreePreparationError("invalid dated-tree checksum row")
        checksums[name] = digest
    inventory = {
        path.name for path in directory.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise TreePreparationError("dated-tree checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(regular(directory / name)) != digest:
            raise TreePreparationError(f"dated-tree checksum mismatch: {name}")
    tree_path = regular(directory / "dated_tree.mean_ma.tre")
    tip_count = validation.get("dated_tree", {}).get("tip_count")
    if not isinstance(tip_count, int) or tip_count < 2:
        raise TreePreparationError("accepted dated-tree validation has an invalid tip count")
    return validation, tree_path


def read_tip_map(path: Path, source_tips: list[str]) -> tuple[list[dict[str, str]], list[str]]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["tree_tip", "biological_species", "include", "rationale"]
        if reader.fieldnames != expected:
            raise TreePreparationError("invalid PGLS tip-map header")
        rows = list(reader)
    if [row["tree_tip"] for row in rows] != list(dict.fromkeys(row["tree_tip"] for row in rows)):
        raise TreePreparationError("duplicate PGLS source-tree tip")
    if set(row["tree_tip"] for row in rows) != set(source_tips):
        raise TreePreparationError("PGLS tip map does not close to the dated tree")
    included: list[str] = []
    for row in rows:
        if row["include"] not in {"true", "false"} or not row["biological_species"] or not row["rationale"]:
            raise TreePreparationError("invalid PGLS tip-map row")
        if row["include"] == "true":
            included.append(row["biological_species"])
    if len(included) < 6 or len(included) != len(set(included)):
        raise TreePreparationError("PGLS included species are insufficient or duplicated")
    return rows, included


def quote_label(label: str) -> str:
    if "'" in label or any(character in label for character in "\t\r\n\x00"):
        raise TreePreparationError(f"unsafe Newick label: {label!r}")
    return f"'{label}'"


def render_newick(clade, root) -> str:
    if clade.is_terminal():
        text = quote_label(str(clade.name))
    else:
        text = "(" + ",".join(render_newick(child, root) for child in clade.clades) + ")"
    if clade is root:
        return text
    return text + f":{float(clade.branch_length):.10f}"


def build_subtree(tree_path: Path, rows: list[dict[str, str]], species_order: list[str]):
    tree = Phylo.read(str(tree_path), "newick")
    terminals = tree.get_terminals()
    source_names = [str(tip.name or "") for tip in terminals]
    selected_rows = [row for row in rows if row["include"] == "true"]
    selected_names = [row["tree_tip"] for row in selected_rows]
    selected_set = set(selected_names)
    source_by_name = {str(tip.name): tip for tip in terminals}
    ancestor = tree.common_ancestor(*(source_by_name[name] for name in selected_names))
    if {str(tip.name) for tip in ancestor.get_terminals()} != selected_set:
        raise TreePreparationError("included PGLS lineages are not one exact dated-tree clade")
    root = copy.deepcopy(ancestor)
    root.branch_length = 0.0
    rename = {row["tree_tip"]: row["biological_species"] for row in selected_rows}
    for tip in root.get_terminals():
        tip.name = rename[str(tip.name)]
    if set(rename) - set(source_names):
        raise TreePreparationError("PGLS tip map names absent source tips")
    internal = [clade for clade in root.find_clades() if not clade.is_terminal()]
    if any(len(clade.clades) != 2 for clade in internal):
        raise TreePreparationError("PGLS subtree is not strictly bifurcating")
    for clade in root.find_clades():
        if clade is root:
            continue
        length = clade.branch_length
        if length is None or not math.isfinite(float(length)) or float(length) < 0:
            raise TreePreparationError("PGLS subtree has invalid branch lengths")
    tip_by_name = {str(tip.name): tip for tip in root.get_terminals()}
    if set(tip_by_name) != set(species_order):
        raise TreePreparationError("PGLS renamed tip set does not close")
    heights = [float(root.distance(tip_by_name[name])) for name in species_order]
    if min(heights) <= 0 or max(heights) - min(heights) > 1e-6:
        raise TreePreparationError("PGLS subtree is not an exact ultrametric Ma tree")
    return root, max(heights) - min(heights), max(heights)


def atomic_publish(output: Path, payloads: dict[str, str]) -> None:
    if os.path.lexists(output):
        raise TreePreparationError(f"refusing to overwrite output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        for name, text in payloads.items():
            (staging / name).write_text(text, encoding="utf-8")
        (staging / "checksums.sha256.tsv").write_text(
            "file\tsha256\n" + "".join(
                f"{name}\t{sha256(staging / name)}\n" for name in sorted(payloads)
            ),
            encoding="utf-8",
        )
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dated-validation-dir", required=True, type=Path)
    parser.add_argument("--tip-map", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        dated_dir = args.dated_validation_dir.expanduser().resolve()
        validation, tree_path = validate_dated_directory(dated_dir)
        source_tree = Phylo.read(str(tree_path), "newick")
        source_tips = [str(tip.name or "") for tip in source_tree.get_terminals()]
        if len(source_tips) != validation["dated_tree"]["tip_count"]:
            raise TreePreparationError("dated-tree tip count differs from its PASS validation")
        rows, species = read_tip_map(args.tip_map, source_tips)
        root, maximum_deviation, height = build_subtree(tree_path, rows, species)
        tree_text = render_newick(root, root) + ":0.0000000000;\n"
        with tempfile.TemporaryDirectory(prefix="pgls-time-tree-check.") as temporary:
            check_path = Path(temporary) / "species_time_tree.nwk"
            check_path.write_text(tree_text, encoding="utf-8")
            check_tree = Phylo.read(str(check_path), "newick")
            if [str(tip.name) for tip in check_tree.get_terminals()] != [
                str(tip.name) for tip in root.get_terminals()
            ]:
                raise TreePreparationError("serialized PGLS tree does not round-trip")
        selected_map = "tree_tip\tbiological_species\n" + "".join(
            f"{row['tree_tip']}\t{row['biological_species']}\n"
            for row in rows if row["include"] == "true"
        )
        output = args.output_dir.expanduser().resolve()
        tree_name = "species_time_tree.nwk"
        map_name = "selected_tip_map.tsv"
        with tempfile.TemporaryDirectory(prefix="pgls-tree-payload.") as temporary:
            temporary_root = Path(temporary)
            temporary_tree = temporary_root / tree_name
            temporary_tree.write_text(tree_text, encoding="utf-8")
            report = {
                "schema_version": "species_time_tree_pass_v1",
                "workflow": "species_time_tree_validation",
                "workflow_version": "1.0.0",
                "status": "PASS",
                "analysis_level": "biological_species",
                "tree": binding(temporary_tree),
                "source_dating_manifest": binding(dated_dir / "validation.json"),
                "biological_species": species,
                "root_semantics": "accepted_biological_species_mrca",
                "branch_length_units": "million_years",
                "checks": {key: True for key in sorted(PASS_CHECKS)},
            }
        atomic_publish(output, {
            tree_name: tree_text,
            map_name: selected_map,
            "species_time_tree_pass.json": json.dumps(
                report, indent=2, sort_keys=True, allow_nan=False
            ) + "\n",
            "tree_metrics.tsv": (
                "metric\tvalue\n"
                f"included_tip_count\t{len(species)}\n"
                f"tree_height_ma\t{height:.12g}\n"
                f"maximum_root_to_tip_deviation_ma\t{maximum_deviation:.12g}\n"
                f"source_root_age_ma\t{validation['dated_tree']['root_age_ma']:.12g}\n"
            ),
        })
        print(f"PASS_SPECIES_PGLS_TIME_TREE\t{len(species)}")
        return 0
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError, TreePreparationError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
