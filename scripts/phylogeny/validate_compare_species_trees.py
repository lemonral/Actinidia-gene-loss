#!/usr/bin/env python3
"""Validate ASTRAL-Pro3 and compare it with a rooted IQ-TREE consensus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Phylo


class ValidationError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    source = path.expanduser().resolve()
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValidationError(f"missing, empty, or symlink file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def expected_terminals(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        labels = {
            row["canonical_tree_label"].strip()
            for row in reader
            if row.get("include_species_tree", "").strip().lower() == "true"
        }
    if not labels or "" in labels:
        raise ValidationError("terminal manifest did not produce valid labels")
    return labels


def terminal_set(tree) -> set[str]:
    labels = {terminal.name for terminal in tree.get_terminals()}
    if None in labels:
        raise ValidationError("tree contains an unnamed terminal")
    return labels


def rooted_clades(tree, all_labels: set[str]) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals(order="preorder"):
        descendants = frozenset(terminal.name for terminal in clade.get_terminals())
        if 1 < len(descendants) < len(all_labels):
            clades.add(descendants)
    return clades


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-tree-state", required=True, type=Path)
    parser.add_argument("--gene-trees", required=True, type=Path)
    parser.add_argument("--astral-tree", required=True, type=Path)
    parser.add_argument("--astral-stderr", required=True, type=Path)
    parser.add_argument("--astral", required=True, type=Path)
    parser.add_argument("--iqtree-consensus", required=True, type=Path)
    parser.add_argument("--terminals", required=True, type=Path)
    parser.add_argument("--root", required=True)
    parser.add_argument("--expected-loci", type=int, default=479)
    parser.add_argument("--expected-version", default="v1.25.3.8")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.output_dir.exists():
            raise ValidationError(f"refusing to overwrite output: {args.output_dir}")
        expected = expected_terminals(args.terminals)
        if args.root not in expected:
            raise ValidationError("declared root is absent from terminal manifest")
        state = json.loads(args.gene_tree_state.read_text(encoding="utf-8"))
        if (
            state.get("status") != "PASS"
            or len(state.get("completed", [])) != args.expected_loci
            or state.get("gene_trees_sha256") != sha256(args.gene_trees)
        ):
            raise ValidationError("gene-tree state/input closure failed")

        tool = args.astral.expanduser().resolve()
        tool_binding = binding(tool)
        probe = subprocess.run([str(tool), "--help"], text=True, capture_output=True)
        banner = probe.stdout + probe.stderr
        if probe.returncode != 0 or args.expected_version not in banner:
            raise ValidationError("ASTRAL-Pro3 version mismatch")
        stderr_text = args.astral_stderr.read_text(encoding="utf-8")
        required = (
            f"Version: {args.expected_version}",
            f"#Genetrees: {args.expected_loci}",
            f"#Species: {len(expected)}",
            "#Duploss: 0",
            "Final Tree:",
            "Score:",
        )
        if any(token not in stderr_text for token in required):
            raise ValidationError("ASTRAL-Pro3 log lacks exact completion evidence")

        astral = Phylo.read(str(args.astral_tree), "newick")
        iqtree = Phylo.read(str(args.iqtree_consensus), "newick")
        if terminal_set(astral) != expected or terminal_set(iqtree) != expected:
            raise ValidationError("species-tree terminal closure failed")
        root_children = astral.root.clades
        if not any(clade.is_terminal() and clade.name == args.root for clade in root_children):
            raise ValidationError("ASTRAL tree is not rooted directly on the declared outgroup")
        iqtree.root_with_outgroup(args.root)
        if not any(clade.is_terminal() and clade.name == args.root for clade in iqtree.root.clades):
            raise ValidationError("IQ-TREE consensus rerooting failed")

        astral_clades = rooted_clades(astral, expected)
        iqtree_clades = rooted_clades(iqtree, expected)
        shared = astral_clades & iqtree_clades
        only_astral = astral_clades - iqtree_clades
        only_iqtree = iqtree_clades - astral_clades

        output = args.output_dir.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            rooted_iqtree = staging / "iqtree_consensus.rooted_vitis.tre"
            Phylo.write(iqtree, str(rooted_iqtree), "newick")
            comparison = staging / "topology_comparison.tsv"
            comparison.write_text(
                "metric\tvalue\n"
                f"astral_internal_clades\t{len(astral_clades)}\n"
                f"iqtree_internal_clades\t{len(iqtree_clades)}\n"
                f"shared_internal_clades\t{len(shared)}\n"
                f"astral_only_internal_clades\t{len(only_astral)}\n"
                f"iqtree_only_internal_clades\t{len(only_iqtree)}\n"
                f"rooted_robinson_foulds_distance\t{len(only_astral) + len(only_iqtree)}\n",
                encoding="utf-8",
            )
            payload = {
                "schema_version": 1,
                "workflow": "astral_pro3_exact_validation_and_topology_comparison",
                "status": "PASS",
                "root": args.root,
                "expected_loci": args.expected_loci,
                "expected_species": len(expected),
                "gene_tree_state": binding(args.gene_tree_state),
                "gene_trees": binding(args.gene_trees),
                "astral": tool_binding,
                "astral_tree": binding(args.astral_tree),
                "astral_stderr": binding(args.astral_stderr),
                "iqtree_consensus": binding(args.iqtree_consensus),
                "terminals": binding(args.terminals),
                "topology": {
                    "astral_internal_clades": len(astral_clades),
                    "iqtree_internal_clades": len(iqtree_clades),
                    "shared_internal_clades": len(shared),
                    "astral_only_internal_clades": len(only_astral),
                    "iqtree_only_internal_clades": len(only_iqtree),
                    "rooted_robinson_foulds_distance": len(only_astral) + len(only_iqtree),
                },
                "checks": {
                    "gene_tree_state_closure": True,
                    "tool_version_closure": True,
                    "clean_completion_evidence": True,
                    "terminal_closure": True,
                    "vitis_root_closure": True,
                },
            }
            validation = staging / "validation.json"
            validation.write_text(
                json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n"
                f"validation.json\t{sha256(validation)}\n"
                f"iqtree_consensus.rooted_vitis.tre\t{sha256(rooted_iqtree)}\n"
                f"topology_comparison.tsv\t{sha256(comparison)}\n",
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(json.dumps(payload["topology"], sort_keys=True))
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, ValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
