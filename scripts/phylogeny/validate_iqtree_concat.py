#!/usr/bin/env python3
"""Validate one completed concatenated IQ-TREE2 analysis fail closed."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from phylo_io import DataError, read_fasta


LEAF = re.compile(r"(?<=[(,])([^():,;\s]+)(?=[:),])")
REQUIRED_SUFFIXES = (
    ".treefile",
    ".contree",
    ".iqtree",
    ".log",
    ".best_scheme.nex",
    ".best_scheme",
    ".splits.nex",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise DataError(f"missing, empty, or symlink input/output: {resolved}")
    return {
        "basename": resolved.name,
        "bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def terminal_labels(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "canonical_tree_label" not in (reader.fieldnames or []):
            raise DataError("terminal manifest lacks canonical_tree_label")
        labels = {
            row["canonical_tree_label"].strip()
            for row in reader
            if row.get("include_species_tree", "").strip().lower() == "true"
        }
    if not labels or "" in labels:
        raise DataError("terminal manifest produced an empty label set")
    return labels


def tree_labels(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8").strip()
    if not text.endswith(";"):
        raise DataError(f"tree is not a complete Newick record: {path.name}")
    return set(LEAF.findall(text))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--partitions", required=True, type=Path)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--iqtree", required=True, type=Path)
    parser.add_argument("--terminals", required=True, type=Path)
    parser.add_argument("--expected-version", default="2.4.0")
    parser.add_argument("--expected-terminal-count", type=int, default=17)
    parser.add_argument("--expected-alignment-length", type=int, default=734085)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    try:
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite validation output: {args.output_dir}")
        inputs = {
            "alignment": binding(args.alignment),
            "partitions": binding(args.partitions),
            "terminals": binding(args.terminals),
        }
        iqtree = args.iqtree.expanduser().resolve()
        tool = binding(iqtree)
        probe = subprocess.run([str(iqtree), "--version"], text=True, capture_output=True)
        version = (probe.stdout + probe.stderr).strip()
        if probe.returncode != 0 or args.expected_version not in version:
            raise DataError(f"IQ-TREE2 version mismatch: {version!r}")

        expected = terminal_labels(args.terminals)
        if len(expected) != args.expected_terminal_count:
            raise DataError("terminal manifest count mismatch")
        alignment = read_fasta(args.alignment)
        if set(alignment) != expected:
            raise DataError("alignment terminal closure failed")
        lengths = {len(sequence) for _, sequence in alignment.values()}
        if lengths != {args.expected_alignment_length}:
            raise DataError(f"alignment length mismatch: {sorted(lengths)}")

        outputs: dict[str, dict[str, object]] = {}
        for suffix in REQUIRED_SUFFIXES:
            path = Path(str(args.prefix) + suffix)
            outputs[suffix.lstrip(".")] = binding(path)
        for suffix in (".treefile", ".contree"):
            path = Path(str(args.prefix) + suffix)
            if tree_labels(path) != expected:
                raise DataError(f"{path.name}: terminal closure failed")
        log_text = Path(str(args.prefix) + ".log").read_text(encoding="utf-8")
        report_text = Path(str(args.prefix) + ".iqtree").read_text(encoding="utf-8")
        if "Date and Time:" not in log_text or "Analysis results written to:" not in log_text:
            raise DataError("IQ-TREE log lacks a clean completion footer")
        if "ultrafast bootstrap" not in report_text.lower():
            raise DataError("IQ-TREE report lacks ultrafast-bootstrap results")
        if "SH-aLRT" not in report_text:
            raise DataError("IQ-TREE report lacks SH-aLRT results")

        output = args.output_dir.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            payload = {
                "schema_version": 1,
                "workflow": "concatenated_iqtree2_exact_validation",
                "status": "PASS",
                "iqtree_version": version,
                "iqtree": tool,
                "inputs": inputs,
                "outputs": outputs,
                "checks": {
                    "alignment_terminal_closure": True,
                    "alignment_length_closure": True,
                    "treefile_terminal_closure": True,
                    "consensus_terminal_closure": True,
                    "clean_completion_footer": True,
                    "ultrafast_bootstrap_present": True,
                    "sh_alrt_present": True,
                },
            }
            validation = staging / "validation.json"
            validation.write_text(
                json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            checksums = staging / "checksums.tsv"
            checksums.write_text(
                "file\tsha256\nvalidation.json\t" + sha256(validation) + "\n",
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS\t{output}")
        return 0
    except (DataError, OSError, subprocess.SubprocessError, UnicodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
