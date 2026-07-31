#!/usr/bin/env python3
"""Trim high-gap codons while preserving reading frame and recording the rule.

The historical script removed nucleotide columns whose gap frequency was
``>= 0.5``.  Because its input was back-translated from a protein alignment,
gaps should occur in complete codons.  This replacement checks that assumption
and trims whole codons, so it cannot quietly turn a codon alignment out of
frame.  It is suitable for a new, versioned run; it does not claim to recover
an undocumented historical command exactly.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_io import DataError, atomic_write_text, read_fasta, report_tsv, safe_output_path, validate_equal_lengths


def trim_one(records: dict[str, tuple[str, str]], path: Path, cutoff: float) -> tuple[str, int, int]:
    length = validate_equal_lengths(records, path)
    if length % 3:
        raise DataError(f"{path}: {length} columns are not divisible by three")
    sequence_count = len(records)
    codon_count = length // 3
    keep_indices: list[int] = []
    for index in range(codon_count):
        start = 3 * index
        gap_characters = 0
        for record_id, (_, sequence) in records.items():
            codon = sequence[start : start + 3]
            if "-" in codon and codon != "---":
                raise DataError(
                    f"{path}: {record_id} has a partial-gap codon at codon {index + 1}: {codon}"
                )
            gap_characters += codon.count("-")
        gap_fraction = gap_characters / (3 * sequence_count)
        if gap_fraction < cutoff:
            keep_indices.append(index)
    if not keep_indices:
        raise DataError(f"{path}: all {codon_count} codons would be removed at cutoff {cutoff}")
    lines: list[str] = []
    for _, (header, sequence) in records.items():
        trimmed = "".join(sequence[index * 3 : index * 3 + 3] for index in keep_indices)
        lines.extend((f">{header}", trimmed))
    return "\n".join(lines) + "\n", codon_count, len(keep_indices)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--glob", default="*.fa")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--gap-cutoff",
        type=float,
        default=0.5,
        help="remove a codon when gap_fraction >= cutoff (legacy rule: 0.5)",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--report", default=None, help="optional TSV report path; default is stdout")
    args = parser.parse_args()
    report: list[dict[str, object]] = []
    try:
        if not 0 < args.gap_cutoff <= 1:
            raise DataError("--gap-cutoff must be >0 and <=1")
        paths = sorted(Path(args.input_dir).glob(args.glob))
        if not paths:
            raise DataError(f"{args.input_dir}: no files match {args.glob!r}")
        output_paths = {path: Path(args.output_dir) / path.name for path in paths}
        for output in output_paths.values():
            safe_output_path(output, args.overwrite)
        rendered: dict[Path, str] = {}
        for path in paths:
            text, input_codons, kept_codons = trim_one(read_fasta(path), path, args.gap_cutoff)
            rendered[output_paths[path]] = text
            report.append(
                {
                    "check": "codon_gap_trim",
                    "severity": "INFO",
                    "status": "PASS",
                    "detail": f"kept {kept_codons}/{input_codons} codons; removed when gap_fraction >= {args.gap_cutoff}",
                    "group": path.name,
                    "output": str(output_paths[path]),
                }
            )
        for output, text in rendered.items():
            atomic_write_text(output, text, overwrite=args.overwrite)
        report.insert(
            0,
            {
                "check": "codon_gap_trim_summary",
                "severity": "INFO",
                "status": "PASS",
                "detail": f"wrote {len(rendered)} frame-preserving trimmed alignments",
            },
        )
        print(report_tsv(report, args.report), end="")
        return 0
    except DataError as error:
        report.append({"check": "codon_gap_trim", "severity": "ERROR", "status": "FAIL", "detail": str(error)})
        print(report_tsv(report, args.report), end="")
        return 2


if __name__ == "__main__":
    sys.exit(main())
