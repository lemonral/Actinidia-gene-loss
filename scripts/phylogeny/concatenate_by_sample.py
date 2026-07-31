#!/usr/bin/env python3
"""Concatenate SCO alignments by explicit gene-to-sample mapping.

The legacy concatenation script assigned the first, second, ... FASTA records
to a hard-coded species list.  That can silently scramble an otherwise valid
matrix if an upstream file changes record order.  This implementation maps each
record through OrthoFinder ``SequenceIDs.txt`` and the reviewed sample manifest,
then writes a supermatrix plus an optional RAxML-NG partition file.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

from phylo_io import (
    DataError,
    atomic_write_text,
    load_samples,
    load_sequence_id_map,
    read_fasta,
    report_tsv,
    resolve_record_samples,
    safe_output_path,
    sample_by_id,
    source_label_to_sample,
    validate_equal_lengths,
)


def partition_name(path: Path) -> str:
    """Return a portable partition identifier based on an orthogroup filename."""

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--glob", default="*.fa")
    parser.add_argument("--sequence-ids", required=True, help="OrthoFinder WorkingDirectory/SequenceIDs.txt")
    parser.add_argument(
        "--species-ids", default=None,
        help="OrthoFinder WorkingDirectory/SpeciesIDs.txt; defaults to the SequenceIDs.txt sibling",
    )
    parser.add_argument("--samples", required=True, help="reviewed samples.tsv")
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-partitions", default=None, help="optional RAxML-NG partition TSV")
    parser.add_argument("--partition-scheme", choices=("gene", "codon"), default="gene")
    parser.add_argument("--expected-groups", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--report", default=None, help="optional TSV report path; default is stdout")
    args = parser.parse_args()

    report: list[dict[str, object]] = []
    try:
        samples = load_samples(args.samples)
        sample_lookup = sample_by_id(samples)
        source_to_sample = source_label_to_sample(samples)
        paths = sorted(Path(args.alignment_dir).glob(args.glob))
        if not paths:
            raise DataError(f"{args.alignment_dir}: no files match {args.glob!r}")
        if args.expected_groups is not None and len(paths) != args.expected_groups:
            raise DataError(f"found {len(paths)} alignments; expected {args.expected_groups}")
        alignments: dict[Path, dict[str, tuple[str, str]]] = {}
        all_ids: set[str] = set()
        for path in paths:
            records = read_fasta(path)
            length = validate_equal_lengths(records, path)
            if length % 3:
                raise DataError(f"{path}: alignment length {length} is not divisible by 3")
            alignments[path] = records
            all_ids.update(records)
        sequence_to_source = load_sequence_id_map(args.sequence_ids, all_ids, args.species_ids)

        output_fasta = safe_output_path(args.output_fasta, args.overwrite)
        output_partitions = None
        if args.output_partitions:
            output_partitions = safe_output_path(args.output_partitions, args.overwrite)
            if output_partitions == output_fasta:
                raise DataError("--output-fasta and --output-partitions must be different files")

        sample_order = [row["sample_id"] for row in samples]
        output_labels = {row["sample_id"]: row["canonical_tree_label"] for row in samples}
        fragments = {sample_id: [] for sample_id in sample_order}
        partitions: list[str] = []
        current_start = 1
        for path in paths:
            records = alignments[path]
            length = validate_equal_lengths(records, path)
            record_to_sample = resolve_record_samples(records, sequence_to_source, source_to_sample)
            sample_to_record: dict[str, str] = {}
            duplicates = [
                sample_id for sample_id, count in Counter(record_to_sample.values()).items() if count > 1
            ]
            if duplicates:
                raise DataError(f"{path}: more than one record maps to {', '.join(sorted(duplicates))}")
            for record_id, sample_id in record_to_sample.items():
                sample_to_record[sample_id] = record_id
            missing = [sample_id for sample_id in sample_order if sample_id not in sample_to_record]
            if missing:
                raise DataError(f"{path}: missing samples: {', '.join(missing)}")
            for sample_id in sample_order:
                fragments[sample_id].append(records[sample_to_record[sample_id]][1])
            current_end = current_start + length - 1
            name = partition_name(path)
            if args.partition_scheme == "gene":
                partitions.append(f"DNA, {name} = {current_start}-{current_end}")
            else:
                for offset in range(3):
                    partitions.append(
                        f"DNA, {name}_pos{offset + 1} = {current_start + offset}-{current_end}\\3"
                    )
            report.append(
                {
                    "check": "concatenate_group",
                    "severity": "INFO",
                    "status": "PASS",
                    "detail": f"{len(records)} mapped records; positions {current_start}-{current_end}",
                    "group": path.name,
                }
            )
            current_start = current_end + 1

        fasta_text = "".join(
            f">{output_labels[sample_id]}\n{''.join(fragments[sample_id])}\n" for sample_id in sample_order
        )
        final_length = len("".join(fragments[sample_order[0]]))
        if any(len("".join(fragments[sample_id])) != final_length for sample_id in sample_order):
            raise DataError("internal error: concatenated sample lengths differ")
        # Every validation has completed before either output is created.
        atomic_write_text(output_fasta, fasta_text, overwrite=args.overwrite)
        if output_partitions is not None:
            atomic_write_text(output_partitions, "\n".join(partitions) + "\n", overwrite=args.overwrite)
        report.insert(
            0,
            {
                "check": "concatenate_summary",
                "severity": "INFO",
                "status": "PASS",
                "detail": f"{len(sample_order)} taxa; {len(paths)} families; {final_length} nt; explicit header mapping",
                "output": str(output_fasta),
            },
        )
        print(report_tsv(report, args.report), end="")
        return 0
    except DataError as error:
        report.append({"check": "concatenate", "severity": "ERROR", "status": "FAIL", "detail": str(error)})
        print(report_tsv(report, args.report), end="")
        return 2


if __name__ == "__main__":
    sys.exit(main())
