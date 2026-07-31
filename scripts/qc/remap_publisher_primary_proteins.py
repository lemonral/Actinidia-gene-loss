#!/usr/bin/env python3
"""Build an exact publisher-primary protein subset with transcript IDs."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.publisher_protein_remap import (
    PublisherProteinRemapError,
    remap_publisher_primary_proteins,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Require complete one-to-one GFF3/publisher protein mapping, subset "
            "to selected primary transcript IDs, and preserve sequences exactly."
        )
    )
    parser.add_argument(
        "--selected-primary-proteins",
        required=True,
        type=Path,
        help="Derived primary protein FASTA; its first-token IDs define the selected set.",
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help=(
            "Complete matched source GFF3 containing transcript/protein accessions; "
            "use full release scope when the publisher protein FASTA is full scope."
        ),
    )
    parser.add_argument(
        "--publisher-proteins",
        required=True,
        type=Path,
        help="Complete publisher protein FASTA from the same release.",
    )
    parser.add_argument(
        "--sample-id", required=True, help="Path-safe assembly-unit identifier."
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New atomic output directory; existing paths are refused.",
    )
    parser.add_argument(
        "--transcript-feature",
        action="append",
        dest="transcript_features",
        help="Accepted transcript feature type; repeat as needed (default: mRNA, transcript).",
    )
    parser.add_argument(
        "--gene-feature",
        action="append",
        dest="gene_features",
        help="Accepted gene-level feature type; repeat as needed (default: gene, pseudogene).",
    )
    parser.add_argument(
        "--gene-as-transcript",
        action="store_true",
        help=(
            "Explicitly accept a gene->CDS graph only when the complete GFF3 has zero "
            "accepted transcript rows. The unchanged gene ID is the self-transcript ID."
        ),
    )
    parser.add_argument(
        "--transcript-id-attribute",
        default="ID",
        help="GFF3 transcript ID attribute (default: ID).",
    )
    parser.add_argument(
        "--transcript-accession-attribute",
        default="Accession",
        help="GFF3 transcript accession attribute (default: Accession).",
    )
    parser.add_argument(
        "--transcript-accession-source",
        choices=("attribute", "transcript_id"),
        default="attribute",
        help=(
            "Read a separate transcript accession attribute (default), or explicitly "
            "use the transcript/self-gene ID itself."
        ),
    )
    parser.add_argument(
        "--cds-parent-attribute",
        default="Parent",
        help="GFF3 CDS parent attribute (default: Parent).",
    )
    parser.add_argument(
        "--protein-accession-attribute",
        default="Protein_Accession",
        help="GFF3 CDS protein accession attribute (default: Protein_Accession).",
    )
    parser.add_argument(
        "--protein-accession-source",
        choices=("attribute", "cds_parent"),
        default="attribute",
        help=(
            "Read the publisher protein accession from a CDS attribute (default), or "
            "explicitly use the sole declared CDS Parent as that accession."
        ),
    )
    parser.add_argument(
        "--publisher-transcript-key",
        default="OriID",
        help="Tab-delimited publisher FASTA transcript header key (default: OriID).",
    )
    parser.add_argument(
        "--publisher-mrna-accession-key",
        default="mRNA",
        help="Tab-delimited publisher FASTA mRNA accession key (default: mRNA).",
    )
    parser.add_argument(
        "--publisher-header-mode",
        choices=("metadata", "first_token"),
        default="metadata",
        help=(
            "Use tab-delimited transcript/mRNA metadata checks (default), or "
            "treat the publisher FASTA first token as the protein accession and "
            "derive transcript mapping only from the GFF3 (NCBI-style bundles)."
        ),
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
    result = remap_publisher_primary_proteins(
        selected_primary_proteins=args.selected_primary_proteins,
        gff_path=args.gff,
        publisher_proteins=args.publisher_proteins,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        transcript_features=args.transcript_features or ("mRNA", "transcript"),
        gene_features=args.gene_features or ("gene", "pseudogene"),
        gene_as_transcript=args.gene_as_transcript,
        transcript_id_attribute=args.transcript_id_attribute,
        transcript_accession_attribute=args.transcript_accession_attribute,
        transcript_accession_source=args.transcript_accession_source,
        cds_parent_attribute=args.cds_parent_attribute,
        protein_accession_attribute=args.protein_accession_attribute,
        protein_accession_source=args.protein_accession_source,
        publisher_transcript_key=args.publisher_transcript_key,
        publisher_mrna_accession_key=args.publisher_mrna_accession_key,
        publisher_header_mode=args.publisher_header_mode,
    )
    print(f"output_dir\t{result.output_dir}")
    print(f"output_proteins\t{result.output_protein_path}")
    print(f"source_publisher_record_count\t{result.source_publisher_record_count}")
    print(f"selected_primary_record_count\t{result.selected_primary_record_count}")
    print(f"excluded_nonprimary_record_count\t{result.excluded_nonprimary_record_count}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PublisherProteinRemapError, OSError, UnicodeError) as error:
        sys.stderr.write(f"error: {error}\n")
        raise SystemExit(2)
