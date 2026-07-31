#!/usr/bin/env python3
"""Fail-closed validation of the sequential TimeTree-bound CAFE5 run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
from io import StringIO
from pathlib import Path

from Bio import Phylo


class ValidationError(RuntimeError):
    pass


MODEL_LAYOUT = {
    "base_poisson": ("Base", ["--poisson"]),
    "gamma3_poisson": ("Gamma", ["--poisson", "--n_gamma_cats", "3"]),
}
CALIBRATION_CLAIM = "TimeTree secondary-calibrated; not fossil-calibrated"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise ValidationError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise ValidationError(f"empty file: {resolved}")
    return resolved


def read_json(path: Path) -> dict[str, object]:
    payload = json.loads(regular(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValidationError(f"JSON object required: {path}")
    return payload


def read_checksum_table(path: Path) -> dict[str, str]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise ValidationError(f"invalid checksum header: {path}")
        rows = list(reader)
    observed: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if (
            not name
            or Path(name).name != name
            or name in observed
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
        ):
            raise ValidationError(f"invalid checksum row: {path}")
        observed[name] = digest
    return observed


def validate_bundle(bundle: Path) -> tuple[dict[str, object], dict[str, str]]:
    if not bundle.is_dir() or bundle.is_symlink():
        raise ValidationError(f"invalid CAFE5 bundle: {bundle}")
    manifest = read_json(bundle / "run_manifest.json")
    if (
        manifest.get("status") != "PASS_PREPARED_CAFE5"
        or manifest.get("workflow") != "cafe5_timetree_secondary_run_bundle"
        or manifest.get("calibration_claim") != CALIBRATION_CLAIM
    ):
        raise ValidationError("CAFE5 run bundle is not the declared PASS bundle")
    checksums = read_checksum_table(bundle / "checksums.tsv")
    inventory = {
        path.name for path in bundle.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise ValidationError("CAFE5 bundle checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(regular(bundle / name)) != digest:
            raise ValidationError(f"CAFE5 bundle checksum mismatch: {name}")
    return manifest, checksums


def read_primary_counts(path: Path, tips: list[str]) -> tuple[list[str], dict[str, list[int]]]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["Desc", "Family ID", *tips]:
            raise ValidationError("primary CAFE5 count header changed")
        rows: dict[str, list[int]] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(tips) + 2 or row[0] != "NA" or not row[1] or row[1] in rows:
                raise ValidationError(f"invalid primary count row: line {line_number}")
            try:
                values = [int(value) for value in row[2:]]
            except ValueError as error:
                raise ValidationError(f"non-integer primary count: line {line_number}") from error
            if any(value < 0 for value in values):
                raise ValidationError(f"negative primary count: line {line_number}")
            rows[row[1]] = values
    return list(rows), rows


def read_analyzed_family_ids(
    count_table: Path, input_family_ids: list[str], console_stdout: Path,
) -> list[str]:
    """Bind CAFE5's declared filtering of families inferred absent at the root."""
    with regular(count_table).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "FamilyID":
            raise ValidationError("missing CAFE5 count-table family header")
        analyzed = [row[0] for row in reader if row]
    if not analyzed or len(analyzed) != len(set(analyzed)):
        raise ValidationError("empty/duplicate analyzed CAFE5 family IDs")
    analyzed_set = set(analyzed)
    if not analyzed_set.issubset(input_family_ids):
        raise ValidationError("CAFE5 analyzed families are not a subset of the frozen input")
    if analyzed != [family_id for family_id in input_family_ids if family_id in analyzed_set]:
        raise ValidationError("CAFE5 analyzed-family order differs from the frozen input order")
    if analyzed != input_family_ids:
        text = regular(console_stdout).read_text(encoding="utf-8")
        matches = re.findall(
            r"^Filtering families not present at the root from:\s*(\d+)\s+to\s+(\d+)\s*$",
            text,
            flags=re.MULTILINE,
        )
        if matches != [(str(len(input_family_ids)), str(len(analyzed)))]:
            raise ValidationError("CAFE5 root-absence filtering log does not close")
    return analyzed


def validate_node_header(header: list[str], tips: list[str]) -> None:
    if len(header) != 2 * len(tips) or header[0] != "FamilyID":
        raise ValidationError("CAFE5 node-table column count is not 2N-1")
    nodes = header[1:]
    if len(nodes) != len(set(nodes)):
        raise ValidationError("duplicate CAFE5 node-table columns")
    for tip in tips:
        if sum(bool(re.fullmatch(re.escape(tip) + r"<\d+>", node)) for node in nodes) != 1:
            raise ValidationError(f"CAFE5 node-table tip closure failed: {tip}")
    if sum(bool(re.fullmatch(r"<\d+>", node)) for node in nodes) != len(tips) - 1:
        raise ValidationError("CAFE5 internal-node column closure failed")


def read_node_table(
    path: Path, tips: list[str], family_ids: list[str], *, nonnegative: bool,
) -> tuple[list[str], dict[str, list[int]]]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValidationError(f"missing CAFE5 node-table header: {path}")
        validate_node_header(header, tips)
        rows: dict[str, list[int]] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header) or not row[0] or row[0] in rows:
                raise ValidationError(f"invalid CAFE5 node row: {path}:{line_number}")
            try:
                values = [int(value) for value in row[1:]]
            except ValueError as error:
                raise ValidationError(f"non-integer CAFE5 node value: {path}:{line_number}") from error
            if nonnegative and any(value < 0 for value in values):
                raise ValidationError(f"negative reconstructed count: {path}:{line_number}")
            rows[row[0]] = values
    if list(rows) != family_ids:
        raise ValidationError(f"CAFE5 family order/closure failed: {path}")
    return header, rows


def read_family_results(path: Path, family_ids: list[str]) -> dict[str, float]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["#FamilyID", "pvalue", "Significant at 0.05"]:
            raise ValidationError(f"unexpected CAFE5 family-result header: {path}")
        rows: dict[str, float] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 3 or not row[0] or row[0] in rows or row[2] not in {"y", "n"}:
                raise ValidationError(f"invalid CAFE5 family-result row: {path}:{line_number}")
            try:
                pvalue = float(row[1])
            except ValueError as error:
                raise ValidationError(f"invalid family p-value: {path}:{line_number}") from error
            if not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
                raise ValidationError(f"non-finite/out-of-range family p-value: {path}:{line_number}")
            if (pvalue < 0.05) != (row[2] == "y"):
                raise ValidationError(f"family significance flag disagrees with p-value: {path}:{line_number}")
            rows[row[0]] = pvalue
    if list(rows) != family_ids:
        raise ValidationError(f"CAFE5 family-result closure failed: {path}")
    return rows


def read_clade_results(path: Path, allowed_nodes: set[str]) -> list[tuple[str, int, int]]:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["#Taxon_ID", "Increase", "Decrease"]:
            raise ValidationError(f"unexpected CAFE5 clade-result header: {path}")
        rows: list[tuple[str, int, int]] = []
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 3 or row[0] not in allowed_nodes or row[0] in seen:
                raise ValidationError(f"invalid CAFE5 clade-result row: {path}:{line_number}")
            try:
                increase, decrease = int(row[1]), int(row[2])
            except ValueError as error:
                raise ValidationError(f"non-integer clade result: {path}:{line_number}") from error
            if increase < 0 or decrease < 0:
                raise ValidationError(f"negative clade result: {path}:{line_number}")
            seen.add(row[0])
            rows.append((row[0], increase, decrease))
    return rows


def validate_family_likelihoods(path: Path, family_ids: set[str]) -> None:
    with regular(path).open(encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines or not lines[0].startswith("#FamilyID\t"):
        raise ValidationError(f"unexpected family-likelihood header: {path}")
    observed: set[str] = set()
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if not fields or fields[0] not in family_ids:
            raise ValidationError(f"unknown likelihood family: {path}:{line_number}")
        observed.add(fields[0])
    if observed != family_ids:
        raise ValidationError(f"family-likelihood closure failed: {path}")


def validate_branch_probabilities(path: Path, node_header: list[str], family_ids: set[str]) -> None:
    with regular(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != node_header:
            raise ValidationError(f"branch-probability header differs from node tables: {path}")
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(node_header) or row[0] not in family_ids or row[0] in seen:
                raise ValidationError(f"invalid branch-probability row: {path}:{line_number}")
            for value in row[1:]:
                if value == "N/A":
                    continue
                try:
                    probability = float(value)
                except ValueError as error:
                    raise ValidationError(f"invalid branch probability: {path}:{line_number}") from error
                if not math.isfinite(probability) or not 0 <= probability <= 1:
                    raise ValidationError(f"out-of-range branch probability: {path}:{line_number}")
            seen.add(row[0])


def validate_asr(path: Path, family_ids: set[str], *, require_all: bool) -> None:
    observed: set[str] = set()
    saw_nexus = saw_begin = saw_end = False
    with regular(path).open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            saw_nexus |= stripped.lower() == "#nexus"
            saw_begin |= stripped.lower() == "begin trees;"
            saw_end |= stripped.lower() == "end;"
            match = re.match(r"TREE\s+(\S+)\s*=", stripped, flags=re.IGNORECASE)
            if match:
                family_id = match.group(1)
                if family_id not in family_ids or family_id in observed:
                    raise ValidationError(f"invalid ASR family ID: {path}")
                observed.add(family_id)
    if not (saw_nexus and saw_begin and saw_end) or not observed:
        raise ValidationError(f"invalid/empty CAFE5 ASR Nexus: {path}")
    if require_all and observed != family_ids:
        raise ValidationError(f"base-model ASR family closure failed: {path}")


def parse_results(path: Path, prefix: str) -> tuple[float, list[float]]:
    text = regular(path).read_text(encoding="utf-8")
    score_match = re.search(
        rf"Model\s+{re.escape(prefix)}\s+(?:Result|Final Likelihood \(-lnL\)):\s*([^\s]+)",
        text,
    )
    lambda_match = re.search(r"^Lambda:\s*(.+)$", text, flags=re.MULTILINE)
    if not score_match or not lambda_match:
        raise ValidationError(f"missing score/lambda in CAFE5 result: {path}")
    try:
        score = float(score_match.group(1))
        lambdas = [float(value) for value in re.split(r"[\s,]+", lambda_match.group(1).strip()) if value]
    except ValueError as error:
        raise ValidationError(f"invalid score/lambda in CAFE5 result: {path}") from error
    if not math.isfinite(score) or not lambdas or any(not math.isfinite(x) or x <= 0 for x in lambdas):
        raise ValidationError(f"non-finite/non-positive score or lambda: {path}")
    return score, lambdas


def tree_signature(tree_text: str) -> tuple[set[str], dict[tuple[str, str], float]]:
    tree = Phylo.read(StringIO(tree_text), "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if not all(tips) or len(tips) != len(set(tips)):
        raise ValidationError("CAFE5 report tree has missing/duplicate tips")
    distances = {
        tuple(sorted((left, right))): tree.distance(left, right)
        for index, left in enumerate(tips) for right in tips[index + 1:]
    }
    return set(tips), distances


def validate_report_tree(path: Path, dated_tree: Path, tips: list[str]) -> list[float]:
    text = regular(path).read_text(encoding="utf-8")
    match = re.search(r"^Tree:(.+)$", text, flags=re.MULTILINE)
    lambda_match = re.search(r"^Lambda:\s*(.+)$", text, flags=re.MULTILINE)
    if not match or not lambda_match:
        raise ValidationError(f"missing tree/lambda in CAFE5 report: {path}")
    report_tips, report_distances = tree_signature(match.group(1).strip())
    dated_tips, dated_distances = tree_signature(regular(dated_tree).read_text(encoding="ascii").strip())
    if report_tips != set(tips) or dated_tips != set(tips) or set(report_distances) != set(dated_distances):
        raise ValidationError(f"CAFE5 report tree tip/pair closure failed: {path}")
    # CAFE5 serializes report-tree branch lengths to about six significant digits.
    for pair, expected in dated_distances.items():
        if not math.isclose(report_distances[pair], expected, rel_tol=5e-6, abs_tol=1e-6):
            raise ValidationError(f"CAFE5 report tree differs from dated tree: {path}")
    try:
        lambdas = [float(value) for value in lambda_match.group(1).split()]
    except ValueError as error:
        raise ValidationError(f"invalid report lambda: {path}") from error
    if not lambdas or any(not math.isfinite(value) or value <= 0 for value in lambdas):
        raise ValidationError(f"non-positive report lambda: {path}")
    return lambdas


def validate_inventory(
    run_dir: Path, row: dict[str, object], bundle: Path, manifest: dict[str, object],
) -> tuple[str, Path]:
    model_id = str(row.get("model_id", ""))
    if model_id not in MODEL_LAYOUT or row.get("returncode") != 0:
        raise ValidationError(f"invalid completed model state: {model_id}")
    prefix, expected_arguments = MODEL_LAYOUT[model_id]
    arguments = next(
        (entry.get("arguments") for entry in manifest["models"] if entry.get("model_id") == model_id),
        None,
    )
    if arguments != expected_arguments:
        raise ValidationError(f"prepared CAFE5 arguments changed: {model_id}")
    results = run_dir / "runs" / model_id / "results"
    if not results.is_dir() or results.is_symlink():
        raise ValidationError(f"invalid result directory: {model_id}")
    expected_command = [
        str(manifest["cafe5"]["path"]), "--infile", str(bundle / "cafe5_primary_lt100.tsv"),
        "--tree", str(bundle / "dated_tree.mean_ma.tre"), "--cores", "1",
        *expected_arguments, "--output_prefix", str(results),
    ]
    if row.get("command") != expected_command:
        raise ValidationError(f"CAFE5 command binding changed: {model_id}")
    for stream in ("stdout", "stderr"):
        stream_path = regular(run_dir / "runs" / model_id / f"console.{stream}", allow_empty=True)
        binding = row.get(stream)
        if not isinstance(binding, dict) or binding.get("bytes") != stream_path.stat().st_size or binding.get("sha256") != sha256(stream_path):
            raise ValidationError(f"CAFE5 {stream} binding changed: {model_id}")
    recorded = row.get("result_files")
    if not isinstance(recorded, list):
        raise ValidationError(f"missing CAFE5 result inventory: {model_id}")
    bindings: dict[str, dict[str, object]] = {}
    for binding in recorded:
        if not isinstance(binding, dict):
            raise ValidationError(f"invalid CAFE5 result binding: {model_id}")
        relative = str(binding.get("relative_path", ""))
        candidate = Path(relative)
        if not relative or candidate.is_absolute() or ".." in candidate.parts or relative in bindings:
            raise ValidationError(f"unsafe/duplicate CAFE5 result path: {model_id}")
        bindings[relative] = binding
    inventory = {
        str(path.relative_to(results)): path for path in results.rglob("*") if path.is_file()
    }
    if set(inventory) != set(bindings):
        raise ValidationError(f"CAFE5 result inventory does not close: {model_id}")
    for relative, path in inventory.items():
        source = regular(path)
        binding = bindings[relative]
        if binding.get("bytes") != source.stat().st_size or binding.get("sha256") != sha256(source):
            raise ValidationError(f"CAFE5 result binding mismatch: {model_id}:{relative}")
    required = {
        f"{prefix}_{suffix}" for suffix in (
            "report.cafe", "results.txt", "family_likelihoods.txt", "asr.tre", "count.tab",
            "change.tab", "family_results.txt", "clade_results.txt", "branch_probabilities.tab",
        )
    }
    if not required.issubset(inventory):
        raise ValidationError(f"missing required CAFE5 results: {model_id}")
    return prefix, results


def atomic_write(output: Path, payloads: dict[str, str]) -> None:
    if output.exists():
        raise ValidationError(f"refusing to overwrite validation output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        for name, text in payloads.items():
            (staging / name).write_text(text, encoding="utf-8")
        names = sorted(payloads)
        (staging / "checksums.tsv").write_text(
            "file\tsha256\n" + "".join(f"{name}\t{sha256(staging / name)}\n" for name in names),
            encoding="utf-8",
        )
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def validate_gamma_initialization_failure(run_dir: Path, state: dict[str, object]) -> dict[str, object]:
    """Bind an unavailable Gamma3 sensitivity without treating it as a successful model."""
    stage = run_dir / "runs" / "gamma3_poisson"
    stdout = regular(stage / "console.stdout")
    stderr = regular(stage / "console.stderr")
    stdout_text = stdout.read_text(encoding="utf-8")
    stderr_text = stderr.read_text(encoding="utf-8")
    if not stdout_text.rstrip().endswith("Failed to initialize any reasonable values"):
        raise ValidationError("Gamma3 failure is not the declared initialization failure")
    required_stderr = (
        "Families with largest size differentials:",
        "You may want to try removing the top few families",
    )
    if any(marker not in stderr_text for marker in required_stderr):
        raise ValidationError("Gamma3 initialization diagnostic is incomplete")
    results = stage / "results"
    if results.exists() and (results.is_symlink() or not results.is_dir() or any(results.iterdir())):
        raise ValidationError("Gamma3 initialization failure unexpectedly produced results")
    expected_missing = str(results / "Gamma_report.cafe")
    if state.get("error") != f"missing or symlink file: {expected_missing}":
        raise ValidationError("Gamma3 controller error does not bind the missing result")
    return {
        "model_id": "gamma3_poisson",
        "status": "UNAVAILABLE_INITIALIZATION_FAILURE",
        "stdout": {"bytes": stdout.stat().st_size, "sha256": sha256(stdout)},
        "stderr": {"bytes": stderr.stat().st_size, "sha256": sha256(stderr)},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--accept-gamma-initialization-failure",
        action="store_true",
        help=(
            "validate the completed Base model while binding an exact Gamma3 initialization "
            "failure as unavailable sensitivity evidence"
        ),
    )
    args = parser.parse_args()
    try:
        bundle = args.bundle.expanduser().resolve()
        run_dir = args.run_dir.expanduser().resolve()
        output = args.output_dir.expanduser().resolve()
        manifest, bundle_checksums = validate_bundle(bundle)
        state = read_json(run_dir / "state.json")
        expected_status = "ERROR" if args.accept_gamma_initialization_failure else "PASS_RUN_COMPLETE"
        if (
            state.get("status") != expected_status
            or state.get("workflow") != "sequential_cafe5_timetree_secondary"
            or state.get("calibration_claim") != CALIBRATION_CLAIM
            or state.get("cores") != 1
        ):
            raise ValidationError("CAFE5 state does not match the requested frozen validation mode")
        if state.get("bundle") != {
            "path": str(bundle),
            "manifest_sha256": bundle_checksums["run_manifest.json"],
            "checksums_sha256": sha256(bundle / "checksums.tsv"),
        }:
            raise ValidationError("CAFE5 state-to-bundle binding changed")
        if state.get("cafe5") != manifest.get("cafe5"):
            raise ValidationError("CAFE5 executable/source metadata changed")
        if state.get("large_family_count_excluded_from_rate_estimation") != manifest.get("large_family_count"):
            raise ValidationError("large-family exclusion count changed")
        if state.get("runtime_outlier_count_excluded_from_rate_estimation", 0) != manifest.get(
            "runtime_outlier_count", 0
        ):
            raise ValidationError("runtime-outlier exclusion count changed")
        completed = state.get("completed")
        expected_completed = ["base_poisson"] if args.accept_gamma_initialization_failure else list(MODEL_LAYOUT)
        if not isinstance(completed, list) or [row.get("model_id") for row in completed] != expected_completed:
            raise ValidationError("CAFE5 models did not finish once in frozen sequential order")
        unavailable_model = (
            validate_gamma_initialization_failure(run_dir, state)
            if args.accept_gamma_initialization_failure else None
        )

        tips = list(manifest.get("tip_order", []))
        if not tips or len(tips) != len(set(tips)) or len(tips) != manifest.get("tip_count"):
            raise ValidationError("invalid frozen CAFE5 tip order")
        family_ids, _ = read_primary_counts(bundle / "cafe5_primary_lt100.tsv", tips)
        if len(family_ids) != manifest.get("primary_family_count"):
            raise ValidationError("primary family count differs from the frozen manifest")
        analyzed_family_ids = read_analyzed_family_ids(
            run_dir / "runs" / "base_poisson" / "results" / "Base_count.tab",
            family_ids,
            run_dir / "runs" / "base_poisson" / "console.stdout",
        )
        family_set = set(analyzed_family_ids)
        model_summaries: list[dict[str, object]] = []
        significant_rows: list[tuple[str, str, float]] = []
        clade_rows: list[tuple[str, str, int, int]] = []
        common_node_header: list[str] | None = None
        for row in completed:
            prefix, results = validate_inventory(run_dir, row, bundle, manifest)
            count_header, _ = read_node_table(
                results / f"{prefix}_count.tab", tips, analyzed_family_ids, nonnegative=True
            )
            change_header, _ = read_node_table(
                results / f"{prefix}_change.tab", tips, analyzed_family_ids, nonnegative=False
            )
            if count_header != change_header or (common_node_header is not None and count_header != common_node_header):
                raise ValidationError("CAFE5 node headers differ between outputs/models")
            common_node_header = count_header
            pvalues = read_family_results(
                results / f"{prefix}_family_results.txt", analyzed_family_ids
            )
            clades = read_clade_results(results / f"{prefix}_clade_results.txt", set(count_header[1:]))
            validate_family_likelihoods(results / f"{prefix}_family_likelihoods.txt", family_set)
            validate_branch_probabilities(
                results / f"{prefix}_branch_probabilities.tab", count_header, family_set
            )
            validate_asr(results / f"{prefix}_asr.tre", family_set, require_all=prefix == "Base")
            score, result_lambdas = parse_results(results / f"{prefix}_results.txt", prefix)
            report_lambdas = validate_report_tree(
                results / f"{prefix}_report.cafe", bundle / "dated_tree.mean_ma.tre", tips
            )
            if len(result_lambdas) != len(report_lambdas) or any(
                not math.isclose(left, right, rel_tol=1e-5, abs_tol=1e-12)
                for left, right in zip(result_lambdas, report_lambdas)
            ):
                raise ValidationError(f"CAFE5 lambda differs between report/results: {prefix}")
            model_id = str(row["model_id"])
            significant = [(family_id, value) for family_id, value in pvalues.items() if value < 0.05]
            significant_rows.extend((model_id, family_id, value) for family_id, value in significant)
            clade_rows.extend((model_id, node, increase, decrease) for node, increase, decrease in clades)
            model_summaries.append({
                "model_id": model_id,
                "role": next(entry["role"] for entry in manifest["models"] if entry["model_id"] == model_id),
                "family_count": len(analyzed_family_ids),
                "significant_family_count_p_lt_0_05": len(significant),
                "score": score,
                "lambda_values": result_lambdas,
                "result_file_count": len(row["result_files"]),
            })

        validation_status = (
            "PASS_CAFE5_BASE_VALIDATED_GAMMA_UNAVAILABLE"
            if args.accept_gamma_initialization_failure else "PASS_CAFE5_VALIDATED"
        )
        validation = {
            "schema_version": 1,
            "workflow": "validated_cafe5_timetree_secondary_models",
            "status": validation_status,
            "calibration_claim": CALIBRATION_CLAIM,
            "bundle": state["bundle"],
            "run_state_sha256": sha256(run_dir / "state.json"),
            "tip_count": len(tips),
            "primary_family_count": len(family_ids),
            "analyzed_family_count_after_root_filter": len(analyzed_family_ids),
            "filtered_not_present_at_root_count": len(family_ids) - len(analyzed_family_ids),
            "large_family_count_excluded_from_rate_estimation": manifest["large_family_count"],
            "runtime_outlier_count_excluded_from_rate_estimation": manifest.get(
                "runtime_outlier_count", 0
            ),
            "models": model_summaries,
        }
        if unavailable_model is not None:
            validation["unavailable_sensitivity"] = unavailable_model
        model_table = "model_id\trole\tfamily_count\tsignificant_family_count_p_lt_0.05\tscore\tlambda_values\tresult_file_count\n" + "".join(
            f"{row['model_id']}\t{row['role']}\t{row['family_count']}\t{row['significant_family_count_p_lt_0_05']}\t{row['score']:.17g}\t{','.join(f'{value:.17g}' for value in row['lambda_values'])}\t{row['result_file_count']}\n"
            for row in model_summaries
        )
        significant_table = "model_id\tfamily_id\tpvalue\n" + "".join(
            f"{model}\t{family}\t{pvalue:.17g}\n" for model, family, pvalue in significant_rows
        )
        clade_table = "model_id\ttaxon_or_node_id\tincrease\tdecrease\n" + "".join(
            f"{model}\t{node}\t{increase}\t{decrease}\n"
            for model, node, increase, decrease in clade_rows
        )
        atomic_write(output, {
            "validation.json": json.dumps(validation, indent=2, sort_keys=True, allow_nan=False) + "\n",
            "model_summary.tsv": model_table,
            "significant_families.tsv": significant_table,
            "clade_summary.tsv": clade_table,
        })
        print(json.dumps(validation, indent=2, sort_keys=True))
        return 0
    except (OSError, UnicodeError, ValueError, KeyError, TypeError, json.JSONDecodeError, ValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
