#!/usr/bin/env python3
"""Fail-closed comparison of derived primary and publisher proteins."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.publisher_protein_qc import (
    PublishedProteinCompatibilityError,
    audit_published_protein_compatibility,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Require exact protein ID-set reconciliation and sequence compatibility "
            "between a derived primary set and its publisher protein FASTA."
        )
    )
    parser.add_argument(
        "--derived-proteins",
        required=True,
        type=Path,
        help="Protein FASTA produced by the primary-annotation standardizer.",
    )
    parser.add_argument(
        "--publisher-proteins",
        required=True,
        type=Path,
        help="Protein FASTA from the same publisher release as the genome and GFF3.",
    )
    parser.add_argument(
        "--sample-id", required=True, help="Path-safe assembly-unit identifier."
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New atomic audit directory; existing paths are refused.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    result = audit_published_protein_compatibility(
        derived_proteins=args.derived_proteins,
        publisher_proteins=args.publisher_proteins,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
    )
    print(f"output_dir\t{result.output_dir}")
    print(f"record_count\t{result.record_count}")
    print(f"exact_record_count\t{result.exact_record_count}")
    print(
        "normalized_exact_record_count\t"
        f"{result.normalized_exact_record_count}"
    )
    print(
        "terminal_stop_normalized_record_count\t"
        f"{result.terminal_stop_normalized_record_count}"
    )
    print(
        "publisher_X_wildcard_record_count\t"
        f"{result.publisher_x_wildcard_record_count}"
    )
    print(
        "publisher_X_wildcard_position_count\t"
        f"{result.publisher_x_wildcard_position_count}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PublishedProteinCompatibilityError, OSError, UnicodeError) as error:
        sys.stderr.write(f"error: {error}\n")
        raise SystemExit(2)
