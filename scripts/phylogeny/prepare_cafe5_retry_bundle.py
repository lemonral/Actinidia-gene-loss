#!/usr/bin/env python3
"""Prepare an audited CAFE5 retry after non-finite initial-rate failure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from copy import deepcopy
from pathlib import Path


class RetryError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise RetryError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise RetryError(f"empty file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path, allow_empty=True)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_checksums(path: Path) -> dict[str, str]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise RetryError("invalid source-bundle checksum header")
        rows = list(reader)
    observed: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if (
            not name or Path(name).name != name or name in observed
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
        ):
            raise RetryError("invalid source-bundle checksum row")
        observed[name] = digest
    return observed


def validate_source_bundle(bundle: Path) -> tuple[dict[str, object], dict[str, str]]:
    if not bundle.is_dir() or bundle.is_symlink():
        raise RetryError(f"invalid source bundle: {bundle}")
    manifest = json.loads(regular(bundle / "run_manifest.json").read_text(encoding="utf-8"))
    if (
        manifest.get("status") != "PASS_PREPARED_CAFE5"
        or manifest.get("workflow") != "cafe5_timetree_secondary_run_bundle"
        or manifest.get("models_run_sequentially") is not True
        or manifest.get("cores_per_model") != 1
        or manifest.get("runtime_outlier_count", 0) != 0
    ):
        raise RetryError("source is not the first exact PASS CAFE5 bundle")
    checksums = read_checksums(bundle / "checksums.tsv")
    inventory = {
        path.name for path in bundle.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise RetryError("source-bundle checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(regular(bundle / name)) != digest:
            raise RetryError(f"source-bundle checksum mismatch: {name}")
    return manifest, checksums


def validate_failed_run(
    failed_run: Path, source_bundle: Path, checksums: dict[str, str],
) -> tuple[dict[str, object], Path, Path, list[tuple[str, int]]]:
    if not failed_run.is_dir() or failed_run.is_symlink():
        raise RetryError(f"invalid failed-run directory: {failed_run}")
    state_path = regular(failed_run / "state.json")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    expected_bundle = {
        "path": str(source_bundle),
        "manifest_sha256": checksums["run_manifest.json"],
        "checksums_sha256": sha256(source_bundle / "checksums.tsv"),
    }
    if (
        state.get("status") != "ERROR"
        or state.get("workflow") != "sequential_cafe5_timetree_secondary"
        or state.get("active_model") != "base_poisson"
        or state.get("completed") != []
        or state.get("bundle") != expected_bundle
        or "Base_report.cafe" not in str(state.get("error", ""))
    ):
        raise RetryError("failed run is not the exact non-finite Base initialization failure")
    stage = failed_run / "runs" / "base_poisson"
    stdout = regular(stage / "console.stdout")
    stderr = regular(stage / "console.stderr")
    stdout_text = stdout.read_text(encoding="utf-8")
    if "Failed to initialize any reasonable values" not in stdout_text:
        raise RetryError("failed CAFE5 stdout lacks the initialization-failure marker")
    scores = re.findall(r"Score \(-lnL\):\s+(\S+)", stdout_text)
    if not scores or any(value.lower() != "inf" for value in scores):
        raise RetryError("failed CAFE5 initialization was not uniformly non-finite")
    warned: list[tuple[str, int]] = []
    for line in stderr.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"(\S+):\s*(\d+)", line.strip())
        if match:
            warned.append((match.group(1), int(match.group(2))))
    if not warned or len({family for family, _ in warned}) != len(warned):
        raise RetryError("failed CAFE5 stderr lacks a unique outlier warning list")
    results = stage / "results"
    if results.exists() and (
        not results.is_dir() or results.is_symlink() or any(results.iterdir())
    ):
        raise RetryError("failed CAFE5 run unexpectedly produced result files")
    return state, stdout, stderr, warned


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-bundle", required=True, type=Path)
    parser.add_argument("--failed-run", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        source = args.source_bundle.expanduser().resolve()
        failed = args.failed_run.expanduser().resolve()
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise RetryError(f"refusing to overwrite output: {output}")
        manifest, checksums = validate_source_bundle(source)
        state, stdout, stderr, warned = validate_failed_run(failed, source, checksums)
        threshold = min(difference for _, difference in warned)
        warned_map = dict(warned)
        with regular(source / "cafe5_primary_lt100.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header != ["Desc", "Family ID", *manifest.get("tip_order", [])]:
                raise RetryError("source primary-count header changed")
            primary: list[list[str]] = []
            outliers: list[list[str]] = []
            observed_differences: dict[str, int] = {}
            for line_number, row in enumerate(reader, start=2):
                if len(row) != len(header) or row[0] != "NA":
                    raise RetryError(f"invalid source family row: line {line_number}")
                try:
                    values = [int(value) for value in row[2:]]
                except ValueError as error:
                    raise RetryError(f"invalid source family count: line {line_number}") from error
                difference = max(values) - min(values)
                observed_differences[row[1]] = difference
                (outliers if difference >= threshold else primary).append(row)
        if not primary or not outliers:
            raise RetryError("runtime-outlier split is empty")
        if len(primary) + len(outliers) != manifest.get("primary_family_count"):
            raise RetryError("retry family split does not close to source primary count")
        for family, difference in warned:
            if observed_differences.get(family) != difference or difference < threshold:
                raise RetryError(f"CAFE5 warning differs from source counts: {family}")
        excluded_ids = {row[1] for row in outliers}
        if not set(warned_map).issubset(excluded_ids):
            raise RetryError("retry threshold does not exclude every warned family")

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            for filename, rows in (
                ("cafe5_primary_lt100.tsv", primary),
                ("cafe5_runtime_outliers.tsv", outliers),
            ):
                with (staging / filename).open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    writer.writerows(rows)
            for filename in ("cafe5_large_ge100.tsv", "dated_tree.mean_ma.tre"):
                shutil.copyfile(regular(source / filename), staging / filename)
            shutil.copyfile(regular(failed / "state.json"), staging / "retry_source_state.json")
            shutil.copyfile(stdout, staging / "retry_source_console.stdout")
            shutil.copyfile(stderr, staging / "retry_source_console.stderr")
            retry_manifest = deepcopy(manifest)
            retry_manifest.update({
                "primary_family_count": len(primary),
                "runtime_outlier_count": len(outliers),
                "runtime_outlier_families_in_rate_estimation": False,
                "rate_estimation_attempt": 2,
                "rate_estimation_family_closure": {
                    "source_primary": manifest["primary_family_count"],
                    "retry_primary": len(primary),
                    "runtime_outliers": len(outliers),
                },
                "retry_after_nonfinite_initialization": {
                    "policy": (
                        "exclude every family whose max-minus-min terminal count is at least "
                        "the minimum difference in CAFE5's deterministic failure warning list"
                    ),
                    "maximum_difference_exclusive": threshold,
                    "warning_family_count": len(warned),
                    "warning_differences": sorted(set(warned_map.values()), reverse=True),
                    "failed_state": binding(staging / "retry_source_state.json"),
                    "failed_stdout": binding(staging / "retry_source_console.stdout"),
                    "failed_stderr": binding(staging / "retry_source_console.stderr"),
                    "source_bundle_manifest_sha256": checksums["run_manifest.json"],
                    "source_bundle_checksums_sha256": sha256(source / "checksums.tsv"),
                    "source_failure_finished_at_utc": state.get("finished_at_utc"),
                },
            })
            (staging / "run_manifest.json").write_text(
                json.dumps(retry_manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            names = sorted(path.name for path in staging.iterdir() if path.is_file())
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n" + "".join(
                    f"{name}\t{sha256(staging / name)}\n" for name in names
                ),
                encoding="utf-8",
            )
            os.replace(staging, output)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(
            f"PASS_PREPARED_CAFE5_RETRY\t{output}\tprimary={len(primary)}\t"
            f"runtime_outliers={len(outliers)}\tdifference_threshold={threshold}"
        )
        return 0
    except (OSError, UnicodeError, ValueError, KeyError, TypeError, json.JSONDecodeError, RetryError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
