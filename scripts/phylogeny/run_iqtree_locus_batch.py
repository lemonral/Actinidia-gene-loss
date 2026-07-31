#!/usr/bin/env python3
"""Run one-thread IQ-TREE2 gene trees sequentially with exact validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from phylo_io import DataError, read_fasta


LEAF = re.compile(r"(?<=[(,])([^():,;\s]+)(?=[:),])")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--glob", default="*.fa")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--iqtree", required=True, type=Path)
    parser.add_argument("--expected-loci", required=True, type=int)
    parser.add_argument("--expected-records", required=True, type=int)
    parser.add_argument("--expected-version", default="2.4.0")
    parser.add_argument("--seed", type=int, default=20260718)
    args = parser.parse_args()

    try:
        iqtree = args.iqtree.resolve()
        if not iqtree.is_file() or not os.access(iqtree, os.X_OK):
            raise DataError(f"IQ-TREE2 is missing or not executable: {iqtree}")
        probe = subprocess.run([str(iqtree), "--version"], text=True, capture_output=True, check=False)
        version = (probe.stdout + probe.stderr).strip()
        if probe.returncode != 0 or args.expected_version not in version:
            raise DataError(f"IQ-TREE2 version mismatch: {version!r}")
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite existing output directory: {args.output_dir}")
        inputs = sorted(args.alignment_dir.glob(args.glob))
        if len(inputs) != args.expected_loci:
            raise DataError(f"found {len(inputs)} loci; expected {args.expected_loci}")
        expected_labels: set[str] | None = None
        input_sha: dict[Path, str] = {}
        for path in inputs:
            records = read_fasta(path)
            if len(records) != args.expected_records:
                raise DataError(f"{path}: found {len(records)} records; expected {args.expected_records}")
            labels = set(records)
            if expected_labels is None:
                expected_labels = labels
            elif labels != expected_labels:
                raise DataError(f"{path}: terminal labels differ from the first locus")
            input_sha[path] = sha256_file(path)
        assert expected_labels is not None

        args.output_dir.mkdir(parents=True)
        runs = args.output_dir / "runs"
        runs.mkdir()
        state: dict[str, object] = {
            "schema_version": 1,
            "status": "running",
            "workflow": "sequential_iqtree2_gene_trees",
            "started_at_utc": now(),
            "expected_loci": len(inputs),
            "completed": [],
            "iqtree_version": version,
            "iqtree_sha256": sha256_file(iqtree),
            "seed": args.seed,
        }
        write_json(args.output_dir / "state.json", state)
        tree_lines: list[str] = []
        completed_rows: list[dict[str, object]] = []
        for index, path in enumerate(inputs):
            locus = path.stem
            run_dir = runs / locus
            run_dir.mkdir()
            prefix = run_dir / locus
            command = [
                str(iqtree), "-s", str(path), "-m", "MFP", "-T", "1",
                "-seed", str(args.seed + index), "--prefix", str(prefix),
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            (run_dir / "console.stdout").write_text(completed.stdout, encoding="utf-8")
            (run_dir / "console.stderr").write_text(completed.stderr, encoding="utf-8")
            if completed.returncode != 0:
                raise DataError(f"{locus}: IQ-TREE2 exited {completed.returncode}")
            treefile = prefix.with_suffix(".treefile")
            iqtree_report = prefix.with_suffix(".iqtree")
            if not treefile.is_file() or not iqtree_report.is_file():
                raise DataError(f"{locus}: required IQ-TREE2 outputs are missing")
            tree = treefile.read_text(encoding="utf-8").strip()
            leaves = set(LEAF.findall(tree))
            if leaves != expected_labels:
                raise DataError(f"{locus}: gene-tree terminal closure failed")
            if sha256_file(path) != input_sha[path]:
                raise DataError(f"{locus}: input changed during IQ-TREE2")
            row = {
                "locus": locus,
                "input_sha256": input_sha[path],
                "tree_sha256": sha256_file(treefile),
                "report_sha256": sha256_file(iqtree_report),
                "seed": args.seed + index,
            }
            completed_rows.append(row)
            tree_lines.append(tree)
            state["completed"] = completed_rows
            state["current_locus"] = locus
            write_json(args.output_dir / "state.json", state)

        trees = args.output_dir / "gene_trees.tre"
        trees.write_text("\n".join(tree_lines) + "\n", encoding="utf-8")
        manifest = args.output_dir / "gene_tree_manifest.tsv"
        columns = ("locus", "input_sha256", "tree_sha256", "report_sha256", "seed")
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in completed_rows:
                handle.write("\t".join(str(row[column]) for column in columns) + "\n")
        if sha256_file(iqtree) != state["iqtree_sha256"]:
            raise DataError("IQ-TREE2 executable changed during the batch")
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": now(),
                "gene_trees_sha256": sha256_file(trees),
                "manifest_sha256": sha256_file(manifest),
            }
        )
        state.pop("current_locus", None)
        write_json(args.output_dir / "state.json", state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return 0
    except (DataError, OSError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
