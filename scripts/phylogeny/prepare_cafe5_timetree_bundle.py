#!/usr/bin/env python3
"""Prepare a CAFE5 bundle bound to the validated TimeTree-secondary dated tree."""

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

from Bio import Phylo


class PreparationError(RuntimeError):
    pass


CALIBRATION_CLAIM = "TimeTree secondary-calibrated; not fossil-calibrated"
LARGE_FAMILY_THRESHOLD = 100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise PreparationError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise PreparationError(f"empty file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_checksum_table(path: Path) -> dict[str, str]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise PreparationError(f"invalid checksum-table header: {path}")
        rows = list(reader)
    observed: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if not name or Path(name).name != name or name in observed:
            raise PreparationError(f"invalid checksum filename: {name!r}")
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise PreparationError(f"invalid checksum digest: {name}")
        observed[name] = digest
    return observed


def validate_count_preparation(directory: Path) -> tuple[dict[str, object], Path]:
    if not directory.is_dir() or directory.is_symlink():
        raise PreparationError(f"invalid CAFE count directory: {directory}")
    summary_path = regular(directory / "preparation.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "PASS_INPUT_PREPARATION_ONLY" or summary.get("workflow") != (
        "topology_bound_cafe5_family_count_preparation"
    ):
        raise PreparationError("CAFE count preparation is not PASS_INPUT_PREPARATION_ONLY")
    checksums = read_checksum_table(directory / "checksums.tsv")
    for name, digest in checksums.items():
        if sha256(regular(directory / name)) != digest:
            raise PreparationError(f"CAFE count checksum mismatch: {name}")
    if checksums.get("preparation.json") != sha256(summary_path):
        raise PreparationError("CAFE count summary is absent from its checksum table")
    matrix = regular(directory / str(summary.get("matrix", {}).get("basename", "")))
    if sha256(matrix) != summary.get("matrix", {}).get("sha256"):
        raise PreparationError("CAFE family matrix differs from its preparation binding")
    return summary, matrix


def validate_dated_tree(validation_path: Path, tree_path: Path) -> tuple[dict[str, object], object]:
    validation = json.loads(regular(validation_path).read_text(encoding="utf-8"))
    if validation.get("status") != "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC":
        raise PreparationError("MCMCTree validation is not PASS_MCMCTREE_VALIDATED_ULTRAMETRIC")
    if validation.get("calibration_claim") != CALIBRATION_CLAIM:
        raise PreparationError("dated-tree calibration claim changed")
    validation_dir = validation_path.expanduser().resolve().parent
    checksums = read_checksum_table(validation_dir / "checksums.tsv")
    inventory = {
        path.name for path in validation_dir.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise PreparationError("MCMCTree validation checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(regular(validation_dir / name)) != digest:
            raise PreparationError(f"MCMCTree validation checksum mismatch: {name}")
    tree_source = regular(tree_path)
    if checksums.get(tree_source.name) != sha256(tree_source):
        raise PreparationError("dated tree is absent from the validation checksum table")
    tree = Phylo.read(str(tree_source), "newick")
    if not tree.rooted:
        tree.rooted = True
    if any(len(clade.clades) != 2 for clade in tree.get_nonterminals()):
        raise PreparationError("CAFE dated tree is not strictly bifurcating")
    if any(clade.branch_length is None or clade.branch_length <= 0 for clade in tree.find_clades() if clade is not tree.root):
        raise PreparationError("CAFE dated tree has missing or non-positive branches")
    distances = [tree.distance(tree.root, tip) for tip in tree.get_terminals()]
    if max(distances) - min(distances) > 1e-6:
        raise PreparationError("CAFE dated tree is not ultrametric")
    return validation, tree


def validate_cafe_install(executable_path: Path) -> dict[str, object]:
    executable = regular(executable_path)
    if not os.access(executable, os.X_OK):
        raise PreparationError(f"CAFE5 is not executable: {executable}")
    root = executable.parent.parent
    source_files = {
        name: regular(root / name)
        for name in ("CHANGELOG.md", "README.md", "LICENSE", "main.cpp", "Makefile", "config.log")
    }
    changelog = source_files["CHANGELOG.md"].read_text(encoding="utf-8")
    config_log = source_files["config.log"].read_text(encoding="utf-8", errors="replace")
    if "## [1.1.0]" not in changelog:
        raise PreparationError("installed CAFE5 source changelog lacks release 1.1.0")
    if "PACKAGE_VERSION='1.1'" not in config_log:
        raise PreparationError("installed CAFE5 build log lacks PACKAGE_VERSION=1.1")
    completed = subprocess.run(
        [str(executable), "--help"], text=True, capture_output=True, check=False, timeout=30
    )
    help_text = completed.stdout + "\n" + completed.stderr
    if completed.returncode != 0 or "Usage: cafe5" not in help_text or "--infile" not in help_text:
        raise PreparationError("CAFE5 help/identity probe failed")
    return {
        "path": str(executable),
        "basename": executable.name,
        "bytes": executable.stat().st_size,
        "sha256": sha256(executable),
        "project_version": "CAFE 5.1.0",
        "embedded_package_version": "1.1",
        "version_evidence": (
            "installed source CHANGELOG [1.1.0] plus config.log PACKAGE_VERSION=1.1; "
            "the executable itself exposes program identity but no version banner"
        ),
        "identity_probe_sha256": hashlib.sha256(help_text.encode("utf-8")).hexdigest(),
        "source_build_bindings": {
            name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for name, path in sorted(source_files.items())
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", required=True, type=Path)
    parser.add_argument("--dated-validation", required=True, type=Path)
    parser.add_argument("--dated-tree", required=True, type=Path)
    parser.add_argument("--cafe5", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise PreparationError(f"refusing to overwrite output: {output}")
        count_summary, matrix = validate_count_preparation(args.counts_dir.expanduser().resolve())
        validation, tree = validate_dated_tree(args.dated_validation, args.dated_tree)
        cafe = validate_cafe_install(args.cafe5)
        tree_tips = [tip.name for tip in tree.get_terminals()]
        if not all(tree_tips) or len(tree_tips) != len(set(tree_tips)):
            raise PreparationError("dated tree has missing or duplicate tips")
        if tree_tips != count_summary.get("terminal_order"):
            raise PreparationError("dated-tree tip order differs from the frozen CAFE count order")

        with matrix.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header != ["Desc", "Family ID", *tree_tips]:
                raise PreparationError("CAFE count header differs from the dated-tree tips")
            primary_rows: list[list[object]] = []
            large_rows: list[list[object]] = []
            seen: set[str] = set()
            maximum = 0
            for line_number, row in enumerate(reader, start=2):
                if len(row) != len(header) or row[0] != "NA" or not row[1] or row[1] in seen:
                    raise PreparationError(f"invalid or duplicate CAFE family row: line {line_number}")
                seen.add(row[1])
                try:
                    counts = [int(value) for value in row[2:]]
                except ValueError as error:
                    raise PreparationError(f"non-integer CAFE count: line {line_number}") from error
                if any(count < 0 or str(count) != value for count, value in zip(counts, row[2:])):
                    raise PreparationError(f"invalid canonical CAFE count: line {line_number}")
                if not any(counts):
                    raise PreparationError(f"all-zero CAFE family: line {line_number}")
                row_maximum = max(counts)
                maximum = max(maximum, row_maximum)
                destination = large_rows if row_maximum >= LARGE_FAMILY_THRESHOLD else primary_rows
                destination.append([row[0], row[1], *counts])
        if len(seen) != count_summary.get("family_count"):
            raise PreparationError("CAFE family row count differs from the frozen preparation")
        if len(large_rows) != count_summary.get("families_with_any_terminal_count_at_least_100"):
            raise PreparationError("large-family count differs from the frozen preparation")
        if not primary_rows:
            raise PreparationError("no CAFE families remain below the documented threshold")

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            for filename, rows in (
                ("cafe5_primary_lt100.tsv", primary_rows),
                ("cafe5_large_ge100.tsv", large_rows),
            ):
                with (staging / filename).open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    writer.writerows(rows)
            shutil.copyfile(regular(args.dated_tree), staging / "dated_tree.mean_ma.tre")
            manifest = {
                "schema_version": 1,
                "workflow": "cafe5_timetree_secondary_run_bundle",
                "status": "PASS_PREPARED_CAFE5",
                "calibration_claim": CALIBRATION_CLAIM,
                "count_preparation": binding(args.counts_dir / "preparation.json"),
                "count_preparation_checksums": binding(args.counts_dir / "checksums.tsv"),
                "source_count_matrix": binding(matrix),
                "dated_validation": binding(args.dated_validation),
                "dated_validation_checksums": binding(args.dated_validation.parent / "checksums.tsv"),
                "dated_tree": binding(staging / "dated_tree.mean_ma.tre"),
                "cafe5": cafe,
                "tip_count": len(tree_tips),
                "tip_order": tree_tips,
                "source_family_count": len(seen),
                "primary_family_count": len(primary_rows),
                "large_family_count": len(large_rows),
                "maximum_observed_family_size": maximum,
                "family_size_policy": {
                    "primary": "maximum terminal count <100",
                    "large_ledger": "maximum terminal count >=100",
                    "rationale": (
                        "installed CAFE5 tutorial states that >=100-copy families can make "
                        "rate estimates non-informative and should be set aside"
                    ),
                    "large_family_followup": "analyse separately only with a fixed primary-model lambda",
                },
                "models": [
                    {
                        "model_id": "base_poisson",
                        "arguments": ["--poisson"],
                        "role": "primary single-rate birth-death model",
                    },
                    {
                        "model_id": "gamma3_poisson",
                        "arguments": ["--poisson", "--n_gamma_cats", "3"],
                        "role": "among-family rate-variation sensitivity",
                    },
                ],
                "cores_per_model": 1,
                "models_run_sequentially": True,
                "large_families_in_rate_estimation": False,
                "dated_tree_validation_summary": {
                    "root_age_ma": validation["dated_tree"]["root_age_ma"],
                    "maximum_root_to_tip_deviation_ma": validation["dated_tree"][
                        "maximum_root_to_tip_deviation_ma"
                    ],
                },
            }
            manifest_path = staging / "run_manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            rows = [
                (path.name, sha256(path)) for path in sorted(staging.iterdir()) if path.is_file()
            ]
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n" + "".join(f"{name}\t{digest}\n" for name, digest in rows),
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS_PREPARED_CAFE5\t{output}")
        return 0
    except (
        OSError, UnicodeError, ValueError, KeyError, json.JSONDecodeError,
        subprocess.SubprocessError, PreparationError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
