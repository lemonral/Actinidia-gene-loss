#!/usr/bin/env python3
"""Run the four frozen minimap2 comparisons used for chromosome numbering.

The runner creates target-to-HY4A, HY4A-to-target, target-to-HY4P, and
HY4P-to-target PAF files sequentially.  It records exact input/output hashes
and never overwrites an existing run directory.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_VERSION = "2.28-r1209"
FIXED_ARGS = ("-x", "asm5", "--secondary=no", "-c", "--cs=long")


class RunError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise RunError(f"{label} is missing or empty: {resolved}")
    return resolved


def fasta_ids(path: Path) -> list[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    identifiers: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split()[0]
                if not identifier:
                    raise RunError(f"Empty FASTA identifier in {path}")
                identifiers.append(identifier)
    if len(identifiers) != 29 or len(set(identifiers)) != 29:
        raise RunError(f"{path} must contain exactly 29 unique chromosome records")
    return identifiers


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> int:
    minimap2 = require_file(args.minimap2, "minimap2 executable")
    if not os.access(minimap2, os.X_OK):
        raise RunError(f"minimap2 is not executable: {minimap2}")
    target = require_file(args.target_genome, "target genome")
    hy4a = require_file(args.hy4a_genome, "HY4A genome")
    hy4p = require_file(args.hy4p_genome, "HY4P genome")
    for path in (target, hy4a, hy4p):
        fasta_ids(path)

    version = subprocess.run(
        [str(minimap2), "--version"], check=True, text=True,
        capture_output=True,
    ).stdout.strip()
    if version != EXPECTED_VERSION:
        raise RunError(f"Expected minimap2 {EXPECTED_VERSION}; found {version}")

    output_dir = args.output_dir.expanduser().resolve()
    try:
        output_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as error:
        raise RunError(f"Output directory already exists: {output_dir}") from error

    inputs = {
        "target": target,
        "hy4a": hy4a,
        "hy4p": hy4p,
    }
    comparisons = (
        ("target_to_hy4a", hy4a, target),
        ("hy4a_to_target", target, hy4a),
        ("target_to_hy4p", hy4p, target),
        ("hy4p_to_target", target, hy4p),
    )
    status_path = output_dir / "status.json"
    provenance = {
        "schema_version": 1,
        "workflow": "bidirectional_chromosome_minimap",
        "unit": args.unit,
        "status": "running",
        "started_at_utc": utc_now(),
        "minimap2_version": version,
        "fixed_argv": ["minimap2", *FIXED_ARGS, "{reference_fasta}", "{query_fasta}"],
        "inputs": {
            key: {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for key, path in inputs.items()
        },
        "comparisons": {},
    }
    atomic_json(status_path, provenance)

    try:
        for name, reference, query in comparisons:
            final_paf = output_dir / f"{name}.paf"
            temporary_paf = output_dir / f".{name}.paf.tmp.{os.getpid()}"
            stderr_path = output_dir / f"{name}.stderr.log"
            command = [str(minimap2), *FIXED_ARGS, str(reference), str(query)]
            started = utc_now()
            with temporary_paf.open("wb") as stdout, stderr_path.open("wb") as stderr:
                completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
            if completed.returncode != 0 or temporary_paf.stat().st_size == 0:
                temporary_paf.unlink(missing_ok=True)
                raise RunError(f"{name} failed with exit code {completed.returncode}")
            os.replace(temporary_paf, final_paf)
            provenance["comparisons"][name] = {
                "reference_role": "hy4a" if reference == hy4a else "hy4p" if reference == hy4p else "target",
                "query_role": "hy4a" if query == hy4a else "hy4p" if query == hy4p else "target",
                "started_at_utc": started,
                "finished_at_utc": utc_now(),
                "exit_code": completed.returncode,
                "paf": {"basename": final_paf.name, "bytes": final_paf.stat().st_size, "sha256": sha256(final_paf)},
                "stderr": {"basename": stderr_path.name, "bytes": stderr_path.stat().st_size, "sha256": sha256(stderr_path)},
            }
            atomic_json(status_path, provenance)
    except BaseException as error:
        provenance["status"] = "failed" if isinstance(error, RunError) else "interrupted"
        provenance["finished_at_utc"] = utc_now()
        provenance["error"] = str(error)
        atomic_json(status_path, provenance)
        raise

    provenance["status"] = "completed"
    provenance["finished_at_utc"] = utc_now()
    atomic_json(status_path, provenance)
    print(f"Completed four chromosome minimap comparisons for {args.unit}")
    return 0


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--unit", required=True)
    p.add_argument("--minimap2", required=True, type=Path)
    p.add_argument("--target-genome", required=True, type=Path)
    p.add_argument("--hy4a-genome", required=True, type=Path)
    p.add_argument("--hy4p-genome", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    return p


def main() -> int:
    try:
        return run(parser().parse_args())
    except (OSError, subprocess.SubprocessError, RunError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
