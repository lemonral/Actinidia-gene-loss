#!/usr/bin/env python3
"""Build topology-bound CAFE5 family counts from validated OrthoFinder memberships."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import Phylo


class CountError(RuntimeError):
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
        raise CountError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def membership_count(cell: str) -> int:
    value = cell.strip()
    if not value:
        return 0
    members = [member.strip() for member in value.split(",")]
    if any(not member for member in members) or len(members) != len(set(members)):
        raise CountError("invalid or duplicate gene IDs within an orthogroup cell")
    return len(members)


def read_terminals(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"terminal_id", "canonical_tree_label", "include_species_tree", "identity_status"}
    if not rows or not required.issubset(rows[0]):
        raise CountError("terminal manifest is empty or missing required columns")
    selected = [row for row in rows if row["include_species_tree"].strip().lower() == "true"]
    if any(row["identity_status"] != "confirmed" for row in selected):
        raise CountError("selected terminal without confirmed identity")
    mapping = {row["canonical_tree_label"]: row["terminal_id"] for row in selected}
    if len(mapping) != len(selected) or len({row["terminal_id"] for row in selected}) != len(selected):
        raise CountError("duplicate selected terminal or canonical tree label")
    return mapping, selected


def validated_orthogroups(
    memberships: Path, validation_path: Path, checksums_path: Path
) -> tuple[dict[str, object], dict[str, str]]:
    validation = json.loads(regular(validation_path).read_text(encoding="utf-8"))
    if validation.get("status") != "PASS" or validation.get("workflow") != "orthofinder3_exact_run_validation":
        raise CountError("OrthoFinder exact validation is not PASS")
    checksums = regular(checksums_path)
    if validation.get("key_output_checksums_sha256") != sha256(checksums):
        raise CountError("key-output checksum table does not match OrthoFinder validation")
    with checksums.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {row["basename"]: row["sha256"] for row in rows}
    source = regular(memberships)
    if expected.get(source.name) != sha256(source):
        raise CountError("Orthogroups membership file is absent from or mismatches validated outputs")
    return validation, expected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--orthofinder-validation", required=True, type=Path)
    parser.add_argument("--orthofinder-checksums", required=True, type=Path)
    parser.add_argument("--terminal-manifest", required=True, type=Path)
    parser.add_argument("--topology-freeze", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise CountError(f"refusing to overwrite output: {output}")
        validation, _ = validated_orthogroups(
            args.orthogroups, args.orthofinder_validation, args.orthofinder_checksums
        )
        label_to_terminal, selected = read_terminals(args.terminal_manifest)

        freeze = json.loads(regular(args.topology_freeze).read_text(encoding="utf-8"))
        topology = regular(args.topology)
        if freeze.get("status") != "PASS" or freeze.get("accepted_topology", {}).get("sha256") != sha256(topology):
            raise CountError("accepted topology is not bound to a PASS topology freeze")
        tree = Phylo.read(str(topology), "newick")
        tree_labels = [tip.name for tip in tree.get_terminals()]
        if set(tree_labels) != set(label_to_terminal) or len(tree_labels) != len(set(tree_labels)):
            raise CountError("tree tips do not close exactly to selected terminal labels")
        ordered_terminal_ids = [label_to_terminal[label] for label in tree_labels]

        source = regular(args.orthogroups)
        family_rows: list[tuple[str, list[int]]] = []
        seen: set[str] = set()
        assigned_by_terminal = {terminal: 0 for terminal in ordered_terminal_ids}
        high_count_100 = 0
        with source.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
                raise CountError("unexpected Orthogroups.tsv header")
            observed_terminals = reader.fieldnames[1:]
            if set(observed_terminals) != set(ordered_terminal_ids):
                raise CountError("Orthogroups.tsv species columns do not close to selected terminals")
            for row in reader:
                family = row["Orthogroup"]
                if not family or family in seen:
                    raise CountError("empty or duplicate orthogroup ID")
                seen.add(family)
                counts = [membership_count(row[terminal]) for terminal in ordered_terminal_ids]
                for terminal, count in zip(ordered_terminal_ids, counts):
                    assigned_by_terminal[terminal] += count
                high_count_100 += int(max(counts) >= 100)
                family_rows.append((family, counts))

        if len(family_rows) != validation.get("orthogroups"):
            raise CountError("orthogroup row count does not match exact validation")
        assigned_total = sum(assigned_by_terminal.values())
        if assigned_total != validation.get("assigned_genes"):
            raise CountError("assigned gene total does not match exact validation")
        if len(selected) != validation.get("species"):
            raise CountError("terminal count does not match exact validation")

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            matrix = staging / "cafe5_family_counts.tsv"
            with matrix.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["Desc", "Family ID", *tree_labels])
                for family, counts in family_rows:
                    writer.writerow(["NA", family, *counts])
            summary = {
                "schema_version": 1,
                "workflow": "topology_bound_cafe5_family_count_preparation",
                "status": "PASS_INPUT_PREPARATION_ONLY",
                "cafe5_execution_allowed": False,
                "execution_gate": "requires matching dated ultrametric tree; current production dating gate is blocked",
                "orthogroups": binding(source),
                "orthofinder_validation": binding(args.orthofinder_validation),
                "orthofinder_checksums": binding(args.orthofinder_checksums),
                "terminal_manifest": binding(args.terminal_manifest),
                "topology_freeze": binding(args.topology_freeze),
                "topology": binding(topology),
                "matrix": binding(matrix),
                "terminal_count": len(tree_labels),
                "terminal_order": tree_labels,
                "family_count": len(family_rows),
                "assigned_gene_count": assigned_total,
                "families_with_any_terminal_count_at_least_100": high_count_100,
                "filtering": "none; raw validated orthogroup counts retained and no arbitrary family-size threshold applied",
            }
            summary_path = staging / "preparation.json"
            summary_path.write_text(
                json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            files = sorted(path for path in staging.iterdir() if path.is_file())
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n" + "".join(f"{path.name}\t{sha256(path)}\n" for path in files),
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS_INPUT_PREPARATION_ONLY\t{output}")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, CountError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
