#!/usr/bin/env python3
"""Build an atomic, audited primary-CDS/protein bundle from FASTA + GFF3.

The program uses one process and explicitly limits numerical-library worker
variables to one.  Production runs should pass ``--require-gffread`` so the
selected IDs, CDS sequences, and protein sequences are independently checked.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.primary_annotation import (
    PrimaryAnnotationError,
    parse_canonical_rule,
    standardize_primary_annotation,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select one validated coding transcript per gene deterministically, "
            "extract CDS/protein sequences, and publish a fail-closed audit bundle."
        )
    )
    parser.add_argument("--genome", required=True, type=Path, help="Chromosome-scope FASTA.")
    parser.add_argument("--gff", required=True, type=Path, help="Matching chromosome-scope GFF3.")
    parser.add_argument("--sample-id", required=True, help="Path-safe assembly-unit identifier.")
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="New directory; existing paths are refused."
    )
    parser.add_argument(
        "--canonical-tag",
        action="append",
        default=[],
        metavar="ATTRIBUTE[=VALUE]",
        help=(
            "Ordered transcript-attribute preference. Repeat for lower-priority rules. "
            "Within a matching rule, the longest valid CDS wins."
        ),
    )
    parser.add_argument(
        "--transcript-feature",
        action="append",
        dest="transcript_features",
        help="Accepted transcript feature type. Repeat as needed (default: mRNA, transcript).",
    )
    parser.add_argument(
        "--gene-feature",
        action="append",
        dest="gene_features",
        help="Accepted gene-level feature type. Repeat as needed (default: gene, pseudogene).",
    )
    parser.add_argument(
        "--missing-phase-policy",
        choices=("fail", "zero"),
        default="fail",
        help=(
            "Reject CDS rows with '.' phase (default), or explicitly treat them as phase 0 and flag them."
        ),
    )
    parser.add_argument(
        "--invalid-coding-gene-policy",
        choices=("fail", "omit"),
        default="fail",
        help=(
            "Fail if a gene has CDS rows but no valid coding isoform (default), or explicitly omit and audit it."
        ),
    )
    parser.add_argument(
        "--gene-as-transcript",
        action="store_true",
        help=(
            "Explicitly permit a gene->CDS publisher graph only when the GFF3 has no "
            "accepted transcript rows. Each declared gene becomes one self-transcript "
            "with its unchanged gene ID; mixed or ambiguous graphs are rejected."
        ),
    )
    parser.add_argument(
        "--gffread",
        default="auto",
        help="gffread executable, 'auto' for PATH discovery (default), or 'none' to skip explicitly.",
    )
    parser.add_argument(
        "--require-gffread",
        action="store_true",
        help="Refuse publication unless gffread is found and every selected sequence agrees.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.require_gffread and args.gffread == "none":
        raise PrimaryAnnotationError("--require-gffread cannot be combined with --gffread none")
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"
    rules = tuple(parse_canonical_rule(value) for value in args.canonical_tag)
    result = standardize_primary_annotation(
        genome_path=args.genome,
        gff_path=args.gff,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        canonical_rules=rules,
        transcript_features=args.transcript_features or ("mRNA", "transcript"),
        gene_features=args.gene_features or ("gene", "pseudogene"),
        missing_phase_policy=args.missing_phase_policy,
        invalid_coding_gene_policy=args.invalid_coding_gene_policy,
        gene_as_transcript=args.gene_as_transcript,
        gffread=args.gffread,
        require_gffread=args.require_gffread,
    )
    print(f"output_dir\t{result.output_dir}")
    print(f"source_gene_count\t{result.source_gene_count}")
    print(f"selected_gene_count\t{result.selected_gene_count}")
    print(f"invalid_coding_gene_count\t{result.invalid_coding_gene_count}")
    print(f"gffread_status\t{result.gffread_status}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PrimaryAnnotationError, OSError, UnicodeError) as error:
        sys.stderr.write(f"error: {error}\n")
        raise SystemExit(2)
