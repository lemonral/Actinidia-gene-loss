#!/usr/bin/env python3
"""Run one reproducible JCVI protein-synteny comparison and summarize both sides.

The workflow reproduces the manuscript-era JCVI parameterization while using
descriptive assembly-unit identifiers, bounded workers, immutable inputs, and
an atomic output directory.  One raw-anchor run yields both reference- and
query-centric gene-depth percentages; these two denominators are always
reported together.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Sequence


SAFE_ALIAS = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")


class JcviRunError(RuntimeError):
    """Raised when a JCVI run cannot satisfy its input or output contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise JcviRunError(f"{label} is missing or empty: {resolved}")
    return resolved


def validate_alias(value: str, label: str) -> str:
    if not SAFE_ALIAS.fullmatch(value):
        raise JcviRunError(
            f"{label} must start with a letter and contain at most 64 ASCII letters, "
            "digits, or underscores"
        )
    return value


def build_commands(
    python_bin: str,
    reference_id: str,
    query_id: str,
    threads: int,
) -> dict[str, list[str]]:
    prefix = f"{reference_id}.{query_id}"
    return {
        "ortholog": [
            python_bin,
            "-m",
            "jcvi.compara.catalog",
            "ortholog",
            reference_id,
            query_id,
            "--dbtype=prot",
            "--align_soft=last",
            f"--cpus={threads}",
            "--cscore=0.7",
            "--tandem_Nmax=10",
            "--dist=20",
            "--min_size=4",
            "--no_strip_names",
            "--no_dotplot",
        ],
        "screen": [
            python_bin,
            "-m",
            "jcvi.compara.synteny",
            "screen",
            f"--qbed={reference_id}.bed",
            f"--sbed={query_id}.bed",
            "--minspan=30",
            "--simple",
            f"{prefix}.anchors",
            f"{prefix}.anchors.simple",
        ],
        "depth": [
            python_bin,
            "-m",
            "jcvi.compara.synteny",
            "depth",
            f"--qbed={reference_id}.bed",
            f"--sbed={query_id}.bed",
            f"--depthfile={prefix}.depth.tsv",
            "--histogram",
            f"{prefix}.anchors",
        ],
    }


def command_text(command: Sequence[str]) -> str:
    return " ".join(shlex.quote(value) for value in command)


def run_command(command: list[str], cwd: Path, log_prefix: Path, environment: dict[str, str]) -> None:
    with log_prefix.with_suffix(".stdout.log").open("w", encoding="utf-8") as stdout, log_prefix.with_suffix(
        ".stderr.log"
    ).open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, cwd=cwd, env=environment, stdout=stdout, stderr=stderr)
    if completed.returncode != 0:
        raise JcviRunError(f"Command failed with exit {completed.returncode}: {command_text(command)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-id", required=True, help="JCVI-safe internal alias")
    parser.add_argument("--query-id", required=True, help="JCVI-safe internal alias")
    parser.add_argument("--reference-display-name", required=True)
    parser.add_argument("--query-display-name", required=True)
    parser.add_argument("--query-accession", required=True)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--reference-bed", required=True, type=Path)
    parser.add_argument("--query-protein", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--minimum-block-size", type=int, default=4)
    parser.add_argument("--allowed-reference-bed-only-ids", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> dict[str, Path]:
    validate_alias(args.reference_id, "--reference-id")
    validate_alias(args.query_id, "--query-id")
    if args.reference_id == args.query_id:
        raise JcviRunError("Reference and query aliases must differ")
    if not 1 <= args.threads <= 10:
        raise JcviRunError("--threads must be between 1 and 10")
    if args.minimum_block_size < 1:
        raise JcviRunError("--minimum-block-size must be positive")
    for value, label in (
        (args.reference_display_name, "--reference-display-name"),
        (args.query_display_name, "--query-display-name"),
        (args.query_accession, "--query-accession"),
    ):
        if not value.strip():
            raise JcviRunError(f"{label} cannot be empty")
    inputs = {
        "reference_protein": require_file(args.reference_protein, "Reference protein FASTA"),
        "reference_bed": require_file(args.reference_bed, "Reference BED"),
        "query_protein": require_file(args.query_protein, "Query protein FASTA"),
        "query_bed": require_file(args.query_bed, "Query BED"),
    }
    if args.allowed_reference_bed_only_ids is not None:
        inputs["reference_bed_only_allowlist"] = require_file(
            args.allowed_reference_bed_only_ids, "Reference BED-only allow-list"
        )
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise JcviRunError(f"Refusing to reuse an existing output directory: {output}")
    return inputs


def main() -> None:
    args = parse_args()
    try:
        inputs = validate_args(args)
        commands = build_commands(args.python_bin, args.reference_id, args.query_id, args.threads)
        plan = {
            "reference_id": args.reference_id,
            "reference_display_name": args.reference_display_name,
            "query_id": args.query_id,
            "query_display_name": args.query_display_name,
            "query_accession": args.query_accession,
            "threads": args.threads,
            "commands": {key: command_text(value) for key, value in commands.items()},
            "inputs": {
                key: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                for key, path in inputs.items()
            },
        }
        if args.validate_only:
            print(json.dumps(plan, indent=2, sort_keys=True))
            return

        if not Path(args.python_bin).is_file() and shutil.which(args.python_bin) is None:
            raise JcviRunError(f"Python executable not found: {args.python_bin}")
        if shutil.which("lastal") is None:
            raise JcviRunError("LAST executable lastal is not on PATH")

        output = args.output_dir.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent))
        work = temporary / "work"
        logs = temporary / "logs"
        work.mkdir()
        logs.mkdir()
        for source, name in (
            (inputs["reference_protein"], f"{args.reference_id}.pep"),
            (inputs["reference_bed"], f"{args.reference_id}.bed"),
            (inputs["query_protein"], f"{args.query_id}.pep"),
            (inputs["query_bed"], f"{args.query_id}.bed"),
        ):
            (work / name).symlink_to(source)

        environment = os.environ.copy()
        for variable in THREAD_VARIABLES:
            environment[variable] = "1"
        for name in ("ortholog", "screen", "depth"):
            run_command(commands[name], work, logs / name, environment)

        prefix = f"{args.reference_id}.{args.query_id}"
        anchors = require_file(work / f"{prefix}.anchors", "Raw JCVI anchors")
        require_file(work / f"{prefix}.anchors.simple", "Screened JCVI anchors")
        depth = require_file(work / f"{prefix}.depth.tsv", "JCVI depth table")
        summarizer = Path(__file__).with_name("summarize_jcvi_depth.py").resolve()
        summary_command = [
            args.python_bin,
            str(summarizer),
            "--sample",
            args.query_id,
            "--display-name",
            args.query_display_name,
            "--accession",
            args.query_accession,
            "--reference-protein",
            str(inputs["reference_protein"]),
            "--reference-bed",
            str(inputs["reference_bed"]),
            "--query-protein",
            str(inputs["query_protein"]),
            "--query-bed",
            str(inputs["query_bed"]),
            "--anchors",
            str(anchors),
            "--depthfile",
            str(depth),
            "--output-tsv",
            str(temporary / "jcvi_bidirectional_coverage.tsv"),
            "--output-json",
            str(temporary / "jcvi_bidirectional_coverage.json"),
            "--minimum-block-size",
            str(args.minimum_block_size),
        ]
        if "reference_bed_only_allowlist" in inputs:
            summary_command.extend(
                ["--allowed-reference-bed-only-ids", str(inputs["reference_bed_only_allowlist"])]
            )
        plan["commands"]["summary"] = command_text(summary_command)
        run_command(summary_command, work, logs / "summary", environment)
        (temporary / "run_manifest.json").write_text(
            json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, output)
        print(f"Completed JCVI comparison: {output}")
    except (OSError, JcviRunError, subprocess.SubprocessError) as exc:
        raise SystemExit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
