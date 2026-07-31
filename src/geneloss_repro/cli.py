"""Command-line interface for the reusable gene-loss workflow modules."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .annotation import extract_annotation
from .blast import classify_tblastx, summarize_classification
from .io_utils import SchemaError, concatenate_tsv
from .master import build_loss_master
from .plotting import plot_loss_summary, plot_spatial_bubble
from .spatial import spatial_summary
from .statistics import ploidy_comparison, subgenome_comparison
from .synorth import call_candidates, normalize_synorth


def _path(value: str) -> Path:
    return Path(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="geneloss",
        description="Portable, auditable utilities for the Actinidia gene-loss workflow.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    concatenate = sub.add_parser("concat-tsv", help="Concatenate TSVs only if their headers are exactly identical.")
    concatenate.add_argument("--input", required=True, type=_path, action="append", help="Repeat once per file, in declared sample order.")
    concatenate.add_argument("--output", required=True, type=_path)

    annotation = sub.add_parser("extract-annotation", help="Extract one primary isoform per gene from genome FASTA + GFF3.")
    annotation.add_argument("--genome", required=True, type=_path)
    annotation.add_argument("--gff", required=True, type=_path)
    annotation.add_argument("--sample-id", required=True)
    annotation.add_argument("--output-dir", required=True, type=_path)

    normalize = sub.add_parser("normalize-synorth", help="Normalize raw SynOrths evidence into a headered TSV.")
    normalize.add_argument("--synorth", required=True, type=_path, help="Raw SynOrths text output retained unchanged elsewhere.")
    normalize.add_argument("--reference-coords", required=True, type=_path)
    normalize.add_argument("--target-sample", required=True)
    normalize.add_argument("--output", required=True, type=_path)
    normalize.add_argument("--reference-side", choices=["auto", "first", "second"], default="auto")

    candidates = sub.add_parser("call-candidates", help="Make putative loss candidates from normalized SynOrths anchors.")
    candidates.add_argument("--reference-coords", required=True, type=_path)
    candidates.add_argument("--synorth", required=True, type=_path, help="Output from normalize-synorth.")
    candidates.add_argument("--output", required=True, type=_path)
    candidates.add_argument("--flank-genes", type=int, default=20)
    candidates.add_argument("--mode", choices=["bracketed", "legacy-neighbor", "all-unmatched"], default="bracketed")
    candidates.add_argument("--min-anchors-each-side", type=int, default=1)

    blast = sub.add_parser("classify-tblastx", help="Classify candidates using tBLASTX evidence with explicit schema validation.")
    blast.add_argument("--candidates", required=True, type=_path)
    blast.add_argument("--blast", required=True, type=_path)
    blast.add_argument("--output", required=True, type=_path)
    blast.add_argument("--schema-output", required=True, type=_path, help="Provenance record for the detected BLAST schema and thresholds.")
    blast.add_argument("--blast-schema", choices=["auto", "blast12", "legacy6-auto", "legacy6-bitscore-evalue", "legacy6-evalue-bitscore"], default="auto")
    blast.add_argument("--min-identity", type=float, default=50.0)
    blast.add_argument("--min-bitscore", type=float, default=50.0)
    blast.add_argument("--max-evalue", type=float, default=1e-5, help="Strict upper bound; hit evalue must be < this value.")
    blast.add_argument("--min-alignment-length", type=int, default=0)
    blast.add_argument("--strictness", choices=["legacy", "synteny-aware"], default="legacy")
    blast.add_argument("--synteny-padding-bp", type=int, default=0)
    blast.add_argument("--uncertain-ids", type=_path, help="One-column IDs in assembly gaps/missing data; they remain uncertain.")
    blast.add_argument("--query-fasta", type=_path, help="Actual FASTA supplied to tBLASTX; candidate IDs absent from it become uncertain, not deleted.")
    blast.add_argument("--compatibility-lists-dir", type=_path, help="Optional old-style *_decayed_genes.txt output directory.")

    summary = sub.add_parser("summarize-loss", help="Summarize classifications against all reference genes.")
    summary.add_argument("--classification", required=True, type=_path)
    summary.add_argument("--reference-coords", required=True, type=_path)
    summary.add_argument("--output", required=True, type=_path)

    master = sub.add_parser("build-loss-master", help="Join classifications and retained anchors into a full downstream-ready master table.")
    master.add_argument("--reference-coords", required=True, type=_path)
    master.add_argument("--classification", required=True, type=_path)
    master.add_argument("--retained-anchors", required=True, type=_path, help="The *.retained_anchors.tsv sidecar from call-candidates.")
    master.add_argument("--output", required=True, type=_path)
    master.add_argument("--sample-metadata", type=_path, help="Optional metadata with target_haplotype/sample_id and ploidy.")
    master.add_argument("--noncandidate-class", choices=["unassessed", "retained_by_synorth"], default="unassessed")
    master.add_argument("--run-id", required=True, help="Immutable identifier for this canonical table.")

    spatial = sub.add_parser("spatial-summary", help="Calculate pseudogene-fragment distribution by chromosome and bin.")
    spatial.add_argument("--classification", required=True, type=_path)
    spatial.add_argument("--target-gff", required=True, type=_path)
    spatial.add_argument("--output-dir", required=True, type=_path)
    spatial.add_argument("--sample-id")
    spatial.add_argument("--chromosome-lengths", type=_path, help="Headered TSV chromosome,length; recommended over deriving max GFF end.")
    spatial.add_argument("--number-of-bins", type=int, default=5)
    spatial.add_argument("--bin-mode", choices=["equal-width", "legacy-nested-midpoint"], default="equal-width")
    spatial.add_argument("--gene-feature", default="gene")
    spatial.add_argument("--loss-class", default="pseudogenized")

    ploidy = sub.add_parser("ploidy-test", help="Mann–Whitney comparison with a declared species or haplotype analysis unit.")
    ploidy.add_argument("--summary", required=True, type=_path)
    ploidy.add_argument("--metadata", required=True, type=_path)
    ploidy.add_argument("--output-dir", required=True, type=_path)
    ploidy.add_argument("--metric", default="assessed_loss_rate")
    ploidy.add_argument("--unit", choices=["species", "haplotype"], default="species")
    ploidy.add_argument("--aggregation", choices=["mean", "weighted"], default="mean")
    ploidy.add_argument("--polyploid-labels", default="tetraploid,hexaploid")
    ploidy.add_argument("--diploid-label", default="diploid")

    subgenome = sub.add_parser("subgenome-test", help="Paired chromosome-level comparison among subgenomes/haplotypes.")
    subgenome.add_argument("--inter-chromosome", required=True, type=_path)
    subgenome.add_argument("--metadata", required=True, type=_path)
    subgenome.add_argument("--output-dir", required=True, type=_path)
    subgenome.add_argument("--metric", default="loss_fragment_per_target_gene")
    subgenome.add_argument("--method", choices=["auto", "paired-t", "wilcoxon", "rm-anova", "friedman"], default="auto")

    loss_plot = sub.add_parser("plot-loss-summary", help="Draw a portable stacked pseudogene/deletion count figure.")
    loss_plot.add_argument("--summary", required=True, type=_path)
    loss_plot.add_argument("--output", required=True, type=_path)
    loss_plot.add_argument("--title", default="")

    spatial_plot = sub.add_parser("plot-spatial-bubble", help="Draw an inter- or intra-chromosome bubble plot.")
    spatial_plot.add_argument("--spatial", required=True, type=_path)
    spatial_plot.add_argument("--output", required=True, type=_path)
    spatial_plot.add_argument("--mode", choices=["inter", "intra"], default="inter")
    spatial_plot.add_argument("--rate-column", default="loss_fragment_per_target_gene")
    return parser


def _print_outputs(outputs: dict[str, Path] | None) -> None:
    if not outputs:
        return
    for label, path in outputs.items():
        print(f"{label}\t{path}")


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "concat-tsv":
            print(f"combined_tsv\t{concatenate_tsv(args.input, args.output)}")
        elif args.command == "extract-annotation":
            _print_outputs(extract_annotation(args.genome, args.gff, args.output_dir, args.sample_id))
        elif args.command == "normalize-synorth":
            _, side = normalize_synorth(args.synorth, args.reference_coords, args.target_sample, args.output, args.reference_side)
            print(f"reference_side\t{side}\nnormalized_synorth\t{args.output}")
        elif args.command == "call-candidates":
            candidates, retained = call_candidates(args.reference_coords, args.synorth, args.output, args.flank_genes, args.mode, args.min_anchors_each_side)
            print(f"putative_candidates\t{len(candidates)}\nretained_anchors\t{len(retained)}\noutput\t{args.output}")
        elif args.command == "classify-tblastx":
            rows = classify_tblastx(
                args.candidates, args.blast, args.output, args.schema_output, args.blast_schema,
                args.min_identity, args.min_bitscore, args.max_evalue, args.min_alignment_length,
                args.strictness, args.synteny_padding_bp, args.uncertain_ids, args.query_fasta, args.compatibility_lists_dir,
            )
            counts = {name: sum(row["classification"] == name for row in rows) for name in ("pseudogenized", "deleted", "uncertain")}
            print("\n".join(f"{name}\t{count}" for name, count in counts.items()))
        elif args.command == "summarize-loss":
            summarize_classification(args.classification, args.reference_coords, args.output)
            print(f"summary\t{args.output}")
        elif args.command == "build-loss-master":
            rows = build_loss_master(
                args.reference_coords, args.classification, args.retained_anchors, args.output,
                args.noncandidate_class, args.sample_metadata, args.run_id,
            )
            eligible = sum(row["rate_eligible"] == "true" for row in rows)
            print(f"master_rows\t{len(rows)}\nrate_eligible_rows\t{eligible}\noutput\t{args.output}")
        elif args.command == "spatial-summary":
            _print_outputs(spatial_summary(
                args.classification, args.target_gff, args.output_dir, args.sample_id, args.chromosome_lengths,
                args.number_of_bins, args.bin_mode, args.gene_feature, args.loss_class,
            ))
        elif args.command == "ploidy-test":
            labels = tuple(item.strip() for item in args.polyploid_labels.split(",") if item.strip())
            _print_outputs(ploidy_comparison(args.summary, args.metadata, args.output_dir, args.metric, args.unit, args.aggregation, labels, args.diploid_label))
        elif args.command == "subgenome-test":
            _print_outputs(subgenome_comparison(args.inter_chromosome, args.metadata, args.output_dir, args.metric, args.method))
        elif args.command == "plot-loss-summary":
            print(f"figure\t{plot_loss_summary(args.summary, args.output, args.title)}")
        elif args.command == "plot-spatial-bubble":
            print(f"figure\t{plot_spatial_bubble(args.spatial, args.output, args.mode, args.rate_column)}")
        else:  # pragma: no cover - argparse enforces known commands
            parser.error(f"unknown command {args.command!r}")
    except (SchemaError, RuntimeError, FileNotFoundError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
