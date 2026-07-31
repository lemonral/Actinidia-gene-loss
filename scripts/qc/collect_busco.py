#!/usr/bin/env python3
"""Collect BUSCO v5 short summaries into one deterministic TSV table.

The collector searches recursively for files named
``short_summary.specific.*.txt``.  It parses the compact BUSCO notation plus
the provenance comments emitted by BUSCO v5.  No third-party Python packages
are required.

Example
-------
python scripts/assembly_qc/collect_busco.py \
  --busco-root results/assembly_qc/busco/genome/runs \
  --output results/assembly_qc/busco/genome/busco_summary.tsv

Output percentages are written without the percent sign.  ``mode`` is the
literal BUSCO-reported mode (for example, ``euk_genome_met`` or ``proteins``),
whereas ``input_path`` is copied from the short-summary provenance line.
Fields absent from an otherwise valid BUSCO v5 summary are left blank.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


NUMBER = r"(?:\d+(?:\.\d*)?|\.\d+)"
NOTATION_RE = re.compile(
    rf"C:\s*(?P<C>{NUMBER})%\s*"
    rf"\[\s*S:\s*(?P<S>{NUMBER})%\s*,\s*D:\s*(?P<D>{NUMBER})%\s*\]\s*,\s*"
    rf"F:\s*(?P<F>{NUMBER})%\s*,\s*M:\s*(?P<M>{NUMBER})%\s*,\s*n:\s*(?P<n>\d+)",
    flags=re.IGNORECASE,
)
COUNT_RE = re.compile(
    r"^\s*(?P<count>\d+)\s+.*\((?P<code>[CSDFM])\)"
    r"(?P<internal_stops>\s+\(of which \d+ contain internal stop codons\))?\s*$"
)

OUTPUT_COLUMNS = (
    "sample",
    "busco_version",
    "dataset",
    "dataset_creation_date",
    "mode",
    "input_path",
    "C_percent",
    "S_percent",
    "D_percent",
    "F_percent",
    "M_percent",
    "n",
    "C_count",
    "S_count",
    "D_count",
    "F_count",
    "M_count",
    "short_summary_path",
)


class BuscoParseError(RuntimeError):
    """Raised when a BUSCO short summary lacks required v5 metrics."""


@dataclass(frozen=True)
class BuscoSummary:
    """Parsed fields from one BUSCO v5 specific short-summary file."""

    sample: str
    busco_version: str
    dataset: str
    dataset_creation_date: str
    mode: str
    input_path: str
    C_percent: str
    S_percent: str
    D_percent: str
    F_percent: str
    M_percent: str
    n: str
    C_count: str
    S_count: str
    D_count: str
    F_count: str
    M_count: str
    short_summary_path: str

    def as_row(self) -> dict[str, str]:
        """Return a DictWriter-compatible row."""
        return {column: getattr(self, column) for column in OUTPUT_COLUMNS}


def comment_value(line: str, marker: str) -> str:
    """Return text following a BUSCO provenance marker, or an empty string."""
    if marker not in line:
        return ""
    return line.split(marker, 1)[1].strip()


def infer_sample(path: Path, dataset: str) -> str:
    """Infer the BUSCO output name from its conventional summary filename."""
    filename = path.name
    prefix = f"short_summary.specific.{dataset}." if dataset else ""
    if prefix and filename.startswith(prefix) and filename.endswith(".txt"):
        sample = filename[len(prefix) : -len(".txt")]
        if sample:
            return sample
    # Runner output is ``runs/<sample>/short_summary...``; this is a safe
    # fallback if an older BUSCO version uses a slightly different filename.
    return path.parent.name


def parse_short_summary(path: Path) -> BuscoSummary:
    """Parse one BUSCO v5 ``short_summary.specific`` text file."""
    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as error:
        raise BuscoParseError(f"Cannot read {path}: {error}") from error

    notation: re.Match[str] | None = None
    version = ""
    dataset = ""
    creation_date = ""
    mode = ""
    input_path = ""
    counts = {code: "" for code in "CSDFM"}

    for line in lines:
        if notation is None:
            notation = NOTATION_RE.search(line)

        value = comment_value(line, "BUSCO version is:")
        if value:
            version = value

        value = comment_value(line, "The lineage dataset is:")
        if value:
            dataset_match = re.match(r"([^\s(]+)", value)
            if dataset_match:
                dataset = dataset_match.group(1)
            date_match = re.search(r"Creation date:\s*([^,)]+)", value, flags=re.IGNORECASE)
            if date_match:
                creation_date = date_match.group(1).strip()

        value = comment_value(line, "BUSCO was run in mode:")
        if value:
            mode = value

        value = comment_value(line, "Summarized benchmarking in BUSCO notation for file")
        if value:
            input_path = value

        count_match = COUNT_RE.match(line)
        if count_match:
            code = count_match.group("code")
            # BUSCO v5 genome-mode summaries may append this exact diagnostic
            # to the Complete (C) row.  Do not accept it on another category.
            if count_match.group("internal_stops") and code != "C":
                continue
            counts[code] = count_match.group("count")

    if notation is None:
        raise BuscoParseError(
            f"{path}: missing complete BUSCO notation C/S/D/F/M/n; "
            "the run may be incomplete or the file may not be a BUSCO v5 short summary"
        )
    if not dataset:
        # The dataset is also encoded in a standard BUSCO filename.  Recover it
        # when possible, but do not guess from arbitrary dots in the sample.
        filename_match = re.match(r"short_summary\.specific\.([^.]+)\..+\.txt$", path.name)
        if filename_match:
            dataset = filename_match.group(1)
    if not dataset:
        raise BuscoParseError(f"{path}: missing lineage dataset provenance")

    groups = notation.groupdict()
    return BuscoSummary(
        sample=infer_sample(path, dataset),
        busco_version=version,
        dataset=dataset,
        dataset_creation_date=creation_date,
        mode=mode,
        input_path=input_path,
        C_percent=groups["C"],
        S_percent=groups["S"],
        D_percent=groups["D"],
        F_percent=groups["F"],
        M_percent=groups["M"],
        n=groups["n"],
        C_count=counts["C"],
        S_count=counts["S"],
        D_count=counts["D"],
        F_count=counts["F"],
        M_count=counts["M"],
        short_summary_path=str(path.resolve()),
    )


def discover_summaries(busco_root: Path) -> list[Path]:
    """Return specific short-summary files in lexical path order."""
    root = Path(busco_root)
    if not root.exists():
        raise BuscoParseError(f"BUSCO root does not exist: {root}")
    if root.is_file():
        if not root.name.startswith("short_summary.specific.") or root.suffix != ".txt":
            raise BuscoParseError(f"Not a BUSCO specific short-summary file: {root}")
        return [root]
    return sorted(root.rglob("short_summary.specific.*.txt"), key=lambda path: str(path))


def collect(paths: Iterable[Path]) -> list[BuscoSummary]:
    """Parse summaries and return rows in deterministic biological order."""
    summaries = [parse_short_summary(path) for path in paths]
    return sorted(
        summaries,
        key=lambda item: (item.sample, item.dataset, item.mode, item.short_summary_path),
    )


def write_tsv(path: Path, summaries: Iterable[BuscoSummary]) -> None:
    """Write the collected table, including a header when there are no rows."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for summary in summaries:
            writer.writerow(summary.as_row())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--busco-root",
        required=True,
        type=Path,
        help="BUSCO batch root, runs directory, or one specific short-summary file",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = discover_summaries(args.busco_root)
        if not paths:
            raise BuscoParseError(
                f"No short_summary.specific.*.txt files found below {args.busco_root}"
            )
        summaries = collect(paths)
        write_tsv(args.output, summaries)
    except BuscoParseError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Collected {len(summaries)} BUSCO summaries into {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
