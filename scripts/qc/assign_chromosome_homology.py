#!/usr/bin/env python3
"""Assign final chromosome labels from four frozen score matrices."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.chromosome_assignment import (
    ChromosomeAssignmentError,
    assign_chromosome_homology,
)
from geneloss_repro.chromosome_provenance import ChromosomeProvenanceError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate exact 29x29 nucleotide/JCVI matrices against HY4A and HY4P, "
            "solve four independent global one-to-one assignments, and publish a "
            "final chromosome map only when every frozen gate passes. This command "
            "does not run alignments."
        )
    )
    parser.add_argument(
        "--nucleotide-hy4a",
        required=True,
        type=Path,
        help="Precomputed long-form nucleotide score matrix against HY4A.",
    )
    parser.add_argument(
        "--jcvi-hy4a",
        required=True,
        type=Path,
        help="Precomputed long-form JCVI score matrix against HY4A.",
    )
    parser.add_argument(
        "--nucleotide-hy4p",
        required=True,
        type=Path,
        help="Precomputed long-form nucleotide score matrix against HY4P.",
    )
    parser.add_argument(
        "--jcvi-hy4p",
        required=True,
        type=Path,
        help="Precomputed long-form JCVI score matrix against HY4P.",
    )
    parser.add_argument(
        "--nucleotide-hy4a-provenance",
        required=True,
        type=Path,
        help="Validated provenance sidecar for --nucleotide-hy4a.",
    )
    parser.add_argument(
        "--jcvi-hy4a-provenance",
        required=True,
        type=Path,
        help="Validated provenance sidecar for --jcvi-hy4a.",
    )
    parser.add_argument(
        "--nucleotide-hy4p-provenance",
        required=True,
        type=Path,
        help="Validated provenance sidecar for --nucleotide-hy4p.",
    )
    parser.add_argument(
        "--jcvi-hy4p-provenance",
        required=True,
        type=Path,
        help="Validated provenance sidecar for --jcvi-hy4p.",
    )
    parser.add_argument(
        "--parameters",
        required=True,
        type=Path,
        help="Project TOML containing the frozen [chromosome_homology] policy.",
    )
    parser.add_argument(
        "--target-asset-registry",
        required=True,
        type=Path,
        help=(
            "Verified target genome/GFF/protein registry for this exact "
            "assembly unit and chromosome scope."
        ),
    )
    parser.add_argument(
        "--reference-asset-registry",
        required=True,
        type=Path,
        help="Frozen HY4A/HY4P genome/GFF/protein/CDS asset registry.",
    )
    parser.add_argument(
        "--reference-chromosome-map-registry",
        required=True,
        type=Path,
        help="Frozen exact HY4A/HY4P reference-ID to canonical-label registry.",
    )
    parser.add_argument(
        "--assembly-unit-id",
        required=True,
        help="Path-safe assembly-unit identifier; not a biological-species replicate.",
    )
    parser.add_argument(
        "--target-scope-id",
        required=True,
        help="Path-safe identifier for the exact target chromosome scope.",
    )
    parser.add_argument(
        "--trusted-repository-commit",
        required=True,
        help=(
            "Full lowercase Git commit ID of the reviewed policy, schemas, and "
            "assignment code used for this run."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New atomic output directory. Existing paths are never overwritten.",
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
    result = assign_chromosome_homology(
        nucleotide_hy4a=args.nucleotide_hy4a,
        jcvi_hy4a=args.jcvi_hy4a,
        nucleotide_hy4p=args.nucleotide_hy4p,
        jcvi_hy4p=args.jcvi_hy4p,
        nucleotide_hy4a_provenance=args.nucleotide_hy4a_provenance,
        jcvi_hy4a_provenance=args.jcvi_hy4a_provenance,
        nucleotide_hy4p_provenance=args.nucleotide_hy4p_provenance,
        jcvi_hy4p_provenance=args.jcvi_hy4p_provenance,
        parameters=args.parameters,
        target_asset_registry=args.target_asset_registry,
        reference_asset_registry=args.reference_asset_registry,
        reference_chromosome_map_registry=args.reference_chromosome_map_registry,
        assembly_unit_id=args.assembly_unit_id,
        target_scope_id=args.target_scope_id,
        trusted_repository_commit=args.trusted_repository_commit,
        output_dir=args.output_dir,
    )
    print(f"output_dir\t{result.output_dir}")
    print(f"status\t{result.status}")
    print(f"publication_gate\t{result.publication_gate}")
    print(f"failure_states\t{';'.join(result.failure_states)}")
    print(f"final_map_row_count\t{result.final_map_row_count}")
    return 0 if result.publication_gate == "PASS" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        ChromosomeAssignmentError,
        ChromosomeProvenanceError,
        OSError,
        UnicodeError,
    ) as error:
        sys.stderr.write(f"error: {error}\n")
        raise SystemExit(2)
