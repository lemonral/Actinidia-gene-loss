#!/usr/bin/env python3
"""Fail closed unless a BUSCO batch is complete, exact, and integrity-valid.

This post-run gate is intended for analysis queues that must not start after a
failed or partial BUSCO prerequisite.  It requires the collected BUSCO TSV to
contain exactly the samples in one manifest, verifies that every corresponding
specific short summary is parseable, applies the same engine-integrity checks
as :mod:`run_busco_batch`, and confirms that each collected row exactly matches
the validated underlying summary.

For genome-mode batches, ``--require-miniprot`` additionally requires actual
Miniprot logs with exactly one normal ``[M::main] Real time:`` completion
marker.  This is stricter than the reusable runner's compatibility behavior,
which permits a genome BUSCO version that uses a different prediction engine.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from collect_busco import OUTPUT_COLUMNS, BuscoParseError, parse_short_summary
from run_busco_batch import BatchInputError, inspect_run_integrity, read_manifest


class PrerequisiteAuditError(RuntimeError):
    """Raised when a BUSCO batch cannot safely release a dependent queue."""


def read_collected_summary(path: Path) -> dict[str, dict[str, str]]:
    """Read an exact-schema collected table and reject duplicate samples."""
    path = path.expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise PrerequisiteAuditError(f"Collected BUSCO summary is missing or empty: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != OUTPUT_COLUMNS:
                raise PrerequisiteAuditError(
                    f"Collected BUSCO summary has an unexpected schema: {path}"
                )
            rows: dict[str, dict[str, str]] = {}
            for line_number, row in enumerate(reader, start=2):
                sample = (row.get("sample") or "").strip()
                if not sample:
                    raise PrerequisiteAuditError(
                        f"Collected BUSCO summary has an empty sample at line {line_number}: {path}"
                    )
                if sample in rows:
                    raise PrerequisiteAuditError(
                        f"Collected BUSCO summary has duplicate sample {sample!r}: {path}"
                    )
                rows[sample] = dict(row)
    except (OSError, UnicodeError, csv.Error) as error:
        raise PrerequisiteAuditError(f"Cannot read collected BUSCO summary {path}: {error}") from error
    return rows


def audit(args: argparse.Namespace) -> int:
    """Validate manifest membership, freshness, summaries, and engine evidence."""
    manifest = args.manifest.expanduser().resolve()
    run_root = args.run_root.expanduser().resolve()
    collected_path = args.collected_summary.expanduser().resolve()

    if args.newer_than_epoch_seconds is not None:
        if args.newer_than_epoch_seconds < 0:
            raise PrerequisiteAuditError("--newer-than-epoch-seconds must be non-negative")
        if not collected_path.exists():
            raise PrerequisiteAuditError(
                f"Collected BUSCO summary was not created after the queue started: {collected_path}"
            )
        observed_mtime = collected_path.stat().st_mtime
        if observed_mtime <= args.newer_than_epoch_seconds:
            raise PrerequisiteAuditError(
                "Collected BUSCO summary was not refreshed by the prerequisite queue: "
                f"mtime {observed_mtime:.6f} <= {args.newer_than_epoch_seconds}"
            )

    manifest_rows = read_manifest(manifest, args.mode, set())
    expected_samples = [sample for sample, _input_path in manifest_rows]
    expected_set = set(expected_samples)
    collected = read_collected_summary(collected_path)
    collected_set = set(collected)
    if collected_set != expected_set:
        missing = sorted(expected_set - collected_set)
        unexpected = sorted(collected_set - expected_set)
        details: list[str] = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        raise PrerequisiteAuditError(
            "Collected BUSCO sample set does not exactly match the manifest: " + "; ".join(details)
        )

    failures: list[str] = []
    passed: list[tuple[str, str]] = []
    for sample in expected_samples:
        run_dir = run_root / sample
        integrity = inspect_run_integrity(run_dir, args.mode)
        if not integrity.complete or integrity.summary_path is None:
            failures.append(f"{sample}: {integrity.message}")
            continue
        if args.require_miniprot and not integrity.message.startswith("Miniprot integrity passed"):
            failures.append(
                f"{sample}: strict Miniprot evidence was required but unavailable; {integrity.message}"
            )
            continue
        try:
            parsed = parse_short_summary(integrity.summary_path)
        except BuscoParseError as error:
            failures.append(f"{sample}: validated summary became unparseable: {error}")
            continue
        if parsed.sample != sample:
            failures.append(
                f"{sample}: underlying summary identifies sample {parsed.sample!r}: {integrity.summary_path}"
            )
            continue
        parsed_row = parsed.as_row()
        mismatched = [
            column
            for column in OUTPUT_COLUMNS
            if collected[sample].get(column, "") != parsed_row[column]
        ]
        if mismatched:
            failures.append(
                f"{sample}: collected row differs from the integrity-valid summary in columns "
                + ",".join(mismatched)
            )
            continue
        passed.append((sample, integrity.message))

    if failures:
        raise PrerequisiteAuditError("BUSCO prerequisite audit failed:\n  " + "\n  ".join(failures))

    for sample, message in passed:
        print(f"PASS\t{sample}\t{message}")
    print(f"Audited {len(passed)} exact BUSCO prerequisite samples successfully")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--collected-summary", required=True, type=Path)
    parser.add_argument("--mode", choices=("genome", "proteins"), required=True)
    parser.add_argument(
        "--require-miniprot",
        action="store_true",
        help="Require strict Miniprot completion evidence for every selected genome run",
    )
    parser.add_argument(
        "--newer-than-epoch-seconds",
        type=int,
        help="Require the collected TSV mtime to be newer than this pre-queue timestamp",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.require_miniprot and args.mode != "genome":
        print("ERROR: --require-miniprot is valid only with --mode genome", file=sys.stderr)
        return 2
    try:
        return audit(args)
    except (PrerequisiteAuditError, BatchInputError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
