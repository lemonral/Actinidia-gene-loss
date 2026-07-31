#!/usr/bin/env python3
"""Freeze an accepted rooted topology after IQ-TREE/ASTRAL agreement passes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import Phylo


class FreezeError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise FreezeError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology-validation", required=True, type=Path)
    parser.add_argument("--astral-rooted", required=True, type=Path)
    parser.add_argument("--iqtree-rooted", required=True, type=Path)
    parser.add_argument("--iqtree-ml-dual-support", required=True, type=Path)
    parser.add_argument("--root", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.output_dir.exists():
            raise FreezeError(f"refusing to overwrite output: {args.output_dir}")
        validation = json.loads(regular(args.topology_validation).read_text(encoding="utf-8"))
        topology = validation.get("topology", {})
        if (
            validation.get("status") != "PASS"
            or topology.get("rooted_robinson_foulds_distance") != 0
            or topology.get("astral_only_internal_clades") != 0
            or topology.get("iqtree_only_internal_clades") != 0
            or validation.get("root") != args.root
        ):
            raise FreezeError("topology validation does not authorize a freeze")
        if validation.get("astral_tree", {}).get("sha256") != sha256(regular(args.astral_rooted)):
            raise FreezeError("ASTRAL tree does not match validation")

        tree = Phylo.read(str(regular(args.astral_rooted)), "newick")
        if not any(clade.is_terminal() and clade.name == args.root for clade in tree.root.clades):
            raise FreezeError("ASTRAL topology is not rooted on the declared outgroup")
        for clade in tree.find_clades():
            clade.branch_length = None
            clade.confidence = None
            if not clade.is_terminal():
                clade.name = None

        output = args.output_dir.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            topology_only = staging / "accepted_topology.rooted_vitis.tre"
            Phylo.write(tree, str(topology_only), "newick", plain=True)
            astral_copy = staging / "accepted_astral.rooted_vitis.tre"
            iqtree_rooted_copy = staging / "accepted_iqtree_consensus.rooted_vitis.tre"
            iqtree_ml_copy = staging / "accepted_iqtree_ml.dual_support.unrooted.tre"
            shutil.copyfile(args.astral_rooted, astral_copy)
            shutil.copyfile(args.iqtree_rooted, iqtree_rooted_copy)
            shutil.copyfile(args.iqtree_ml_dual_support, iqtree_ml_copy)
            manifest = {
                "schema_version": 1,
                "workflow": "selected_lineage_species_topology_freeze",
                "status": "PASS",
                "root": args.root,
                "decision": "IQ-TREE concatenation and ASTRAL-Pro3 have identical rooted topology",
                "rooted_robinson_foulds_distance": 0,
                "topology_validation": binding(args.topology_validation),
                "accepted_topology": binding(topology_only),
                "astral_rooted": binding(astral_copy),
                "iqtree_consensus_rooted": binding(iqtree_rooted_copy),
                "iqtree_ml_dual_support_unrooted": binding(iqtree_ml_copy),
            }
            manifest_path = staging / "topology_freeze.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            files = sorted(path for path in staging.iterdir() if path.is_file())
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n"
                + "".join(f"{path.name}\t{sha256(path)}\n" for path in files),
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS\t{output}")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, FreezeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
