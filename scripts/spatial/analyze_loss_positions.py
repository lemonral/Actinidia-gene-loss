#!/usr/bin/env python3
"""Analyze positive gene-loss fragment positions across assembly units.

Inputs are deliberately separate and explicit:

* a final call table (complete matrix or positive-only table),
* exactly one selected target feature coordinate per analyzed positive call,
* an assembly-unit manifest with one genome FASTA and GFF per unit, and
* optionally, independently supported centromere intervals.

The primary outputs use mutually exclusive equal-width chromosome bins and a
0-at-end/1-at-center normalized distance.  ``--legacy-reproduction`` writes
the manuscript-era overlapping nested-midpoint intervals as a clearly labelled
sensitivity artifact; those intervals are not valid categories for an
inferential test.

The script has no species-specific sample limit.  Give *A. deliciosa* A-F and
*A. zhejiangensis* A/B separate manifest rows and keep their suffixes in
``haplotype_or_subgenome`` so every assembly unit remains visible in output.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.io_utils import SchemaError
from geneloss_repro.spatial import analyze_loss_positions


def _optional_column(value: str) -> str | None:
    normalized = value.strip()
    return normalized or None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--positive-calls", required=True, type=Path)
    parser.add_argument("--feature-coordinates", required=True, type=Path)
    parser.add_argument("--assembly-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--analysis-label",
        required=True,
        help="Reader-facing scope label, e.g. primary_nonshared_pseudogenized.",
    )
    parser.add_argument(
        "--positive-classes",
        default="pseudogenized",
        help=(
            "Comma-separated final classifications to analyze. Deleted calls normally lack a "
            "target coordinate and should not be added unless an explicit coordinate exists."
        ),
    )
    parser.add_argument("--number-of-bins", type=int, default=5)
    parser.add_argument(
        "--centromeres",
        type=Path,
        help=(
            "Optional TSV: assembly_unit_id, chromosome, centromere_start, centromere_end, "
            "evidence_source. Intervals must be independently supplied."
        ),
    )
    parser.add_argument(
        "--require-complete-centromeres",
        action="store_true",
        help="Fail unless every analyzed GFF gene-bearing sequence has a centromere interval.",
    )
    parser.add_argument(
        "--legacy-reproduction",
        action="store_true",
        help="Also write labelled overlapping manuscript-era nested midpoint intervals.",
    )
    parser.add_argument("--gene-feature", default="gene")

    call = parser.add_argument_group("positive-call column mapping")
    call.add_argument("--call-unit-column", default="target_haplotype")
    call.add_argument("--call-gene-column", default="reference_gene_id")
    call.add_argument("--call-classification-column", default="classification")

    coordinates = parser.add_argument_group("feature-coordinate column mapping")
    coordinates.add_argument("--coordinate-unit-column", default="target_haplotype")
    coordinates.add_argument("--coordinate-gene-column", default="reference_gene_id")
    coordinates.add_argument("--coordinate-chromosome-column", default="target_chromosome")
    coordinates.add_argument("--coordinate-start-column", default="target_start")
    coordinates.add_argument("--coordinate-end-column", default="target_end")
    coordinates.add_argument(
        "--coordinate-classification-column",
        default="classification",
        help=(
            "Set to an empty string only when the coordinate table contains positive calls "
            "exclusively. Otherwise non-positive coordinate rows are filtered by this column."
        ),
    )

    manifest = parser.add_argument_group("assembly-manifest column mapping")
    manifest.add_argument("--manifest-unit-column", default="assembly_unit_id")
    manifest.add_argument("--manifest-species-column", default="biological_species")
    manifest.add_argument("--manifest-haplotype-column", default="haplotype_or_subgenome")
    manifest.add_argument("--manifest-scope-column", default="assembly_scope")
    manifest.add_argument("--manifest-genome-column", default="genome")
    manifest.add_argument("--manifest-gff-column", default="gff")
    manifest.add_argument(
        "--manifest-include-column",
        default="",
        help=(
            "Optional strict true/false column. When set, selected manifest units must equal "
            "the positive-call assembly-unit scope exactly."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    classes = [value.strip() for value in args.positive_classes.split(",") if value.strip()]
    try:
        outputs = analyze_loss_positions(
            args.positive_calls,
            args.feature_coordinates,
            args.assembly_manifest,
            args.output_dir,
            analysis_label=args.analysis_label,
            positive_classes=classes,
            number_of_bins=args.number_of_bins,
            centromeres_path=args.centromeres,
            require_complete_centromeres=args.require_complete_centromeres,
            legacy_reproduction=args.legacy_reproduction,
            gene_feature=args.gene_feature,
            call_unit_column=args.call_unit_column,
            call_gene_column=args.call_gene_column,
            call_classification_column=args.call_classification_column,
            coordinate_unit_column=args.coordinate_unit_column,
            coordinate_gene_column=args.coordinate_gene_column,
            coordinate_chromosome_column=args.coordinate_chromosome_column,
            coordinate_start_column=args.coordinate_start_column,
            coordinate_end_column=args.coordinate_end_column,
            coordinate_classification_column=_optional_column(
                args.coordinate_classification_column
            ),
            manifest_unit_column=args.manifest_unit_column,
            manifest_species_column=args.manifest_species_column,
            manifest_haplotype_column=args.manifest_haplotype_column,
            manifest_scope_column=args.manifest_scope_column,
            manifest_genome_column=args.manifest_genome_column,
            manifest_gff_column=args.manifest_gff_column,
            manifest_include_column=_optional_column(args.manifest_include_column),
        )
    except (OSError, SchemaError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for label, path in outputs.items():
        print(f"{label}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
