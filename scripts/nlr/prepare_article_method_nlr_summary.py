#!/usr/bin/env python3
"""Build a 23-unit NLR summary with classified non-shared loss evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)
RESOLVED = {"retained", "decayed", "deleted"}
POSITIVE = {"decayed", "deleted"}
LOSS_TYPE_ORDER = (
    "no_qualifying_translated_hit",
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
    "truncation_or_partial_alignment_candidate",
    "residual_sequence_mechanism_unresolved",
)
TRUNCATION_CAUSES = {
    "n_terminal_alignment_truncation_candidate",
    "c_terminal_alignment_truncation_candidate",
    "both_terminal_alignment_truncation_candidate",
    "partial_local_alignment_other_candidate",
}
UNRESOLVED_CAUSES = {
    "frameshift_or_stop_below_strict_quality_gate",
    "local_sequence_no_explicit_coding_disruption",
    "genomewide_tblastx_hit_noncallable_local_locus",
    "genomewide_tblastx_hit_without_local_miniprot_support",
}


class NlrSummaryError(ValueError):
    """Raised when the frozen NLR and loss-classification inputs do not reconcile."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article-matrix", required=True, type=Path)
    parser.add_argument("--shared-genes", required=True, type=Path)
    parser.add_argument("--reference-nlr-universe", required=True, type=Path)
    parser.add_argument("--repertoire-counts", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-species", type=int, default=13)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-reference-nlrs", type=int, default=214)
    parser.add_argument("--expected-nonshared-reference-nlrs", type=int, default=76)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise NlrSummaryError(f"missing or empty input: {path.name}")
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    if not fields or not rows or len(fields) != len(set(fields)):
        raise NlrSummaryError(f"invalid TSV input: {path.name}")
    return rows, fields


def require(fields: list[str], needed: set[str], label: str) -> None:
    missing = sorted(needed - set(fields))
    if missing:
        raise NlrSummaryError(f"{label} missing columns: {', '.join(missing)}")


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                writer.writerows(rows)


def loss_type(refined_cause: str) -> str:
    if refined_cause == "no_qualifying_genomewide_tblastx_hit":
        return "no_qualifying_translated_hit"
    if refined_cause == "frameshift_supported":
        return "frameshift_supported"
    if refined_cause == "stop_supported":
        return "inframe_stop_supported"
    if refined_cause == "frameshift_and_stop_supported":
        return "frameshift_and_stop_supported"
    if refined_cause in TRUNCATION_CAUSES:
        return "truncation_or_partial_alignment_candidate"
    if refined_cause in UNRESOLVED_CAUSES:
        return "residual_sequence_mechanism_unresolved"
    raise NlrSummaryError(f"unmapped positive-loss cause: {refined_cause!r}")


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise NlrSummaryError(f"output directory already exists: {args.output_dir}")

    metadata_rows, metadata_fields = read_tsv(args.unit_metadata)
    require(
        metadata_fields,
        {
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "assembly_scope",
            "include",
        },
        args.unit_metadata.name,
    )
    metadata_rows = [
        row for row in metadata_rows if row["include"].strip().lower() == "true"
    ]
    units = [row["assembly_unit_id"] for row in metadata_rows]
    if len(units) != args.expected_units or len(set(units)) != args.expected_units:
        raise NlrSummaryError("unit metadata is not the expected unique cohort")
    if len({row["biological_species"] for row in metadata_rows}) != args.expected_species:
        raise NlrSummaryError("unit metadata species count changed")

    universe_rows, universe_fields = read_tsv(args.reference_nlr_universe)
    require(
        universe_fields,
        {"reference_nlr_id"},
        args.reference_nlr_universe.name,
    )
    reference_nlrs = {row["reference_nlr_id"] for row in universe_rows}
    if (
        len(universe_rows) != args.expected_reference_nlrs
        or len(reference_nlrs) != args.expected_reference_nlrs
        or "" in reference_nlrs
    ):
        raise NlrSummaryError("reference NLR universe changed")

    shared_rows, shared_fields = read_tsv(args.shared_genes)
    require(shared_fields, {"reference_gene_id"}, args.shared_genes.name)
    shared_genes = {row["reference_gene_id"] for row in shared_rows}
    shared_nlrs = reference_nlrs & shared_genes
    nonshared_nlrs = reference_nlrs - shared_genes
    if len(nonshared_nlrs) != args.expected_nonshared_reference_nlrs:
        raise NlrSummaryError(
            "article-method non-shared reference NLR count changed: "
            f"{len(nonshared_nlrs)}"
        )

    repertoire_rows, repertoire_fields = read_tsv(args.repertoire_counts)
    require(
        repertoire_fields,
        {"assembly_unit_id", "total_nlr_count"},
        args.repertoire_counts.name,
    )
    repertoire: dict[str, dict[str, str]] = {}
    for row in repertoire_rows:
        unit = row["assembly_unit_id"]
        if unit in repertoire:
            raise NlrSummaryError(f"duplicate repertoire row: {unit}")
        repertoire[unit] = row
    if set(repertoire) != set(units):
        raise NlrSummaryError("repertoire counts do not match the exact unit cohort")

    classifications: dict[tuple[str, str], dict[str, str]] = {}
    matrix_rows = 0
    with open_text(args.article_matrix) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        require(
            fields,
            {
                "reference_gene_id",
                "assembly_unit_id",
                "manuscript_classification",
                "refined_decayed_cause",
                "refined_cause_evidence_level",
            },
            args.article_matrix.name,
        )
        for line_number, row in enumerate(reader, 2):
            matrix_rows += 1
            gene = row["reference_gene_id"]
            if gene not in reference_nlrs:
                continue
            unit = row["assembly_unit_id"]
            key = (unit, gene)
            if unit not in repertoire or key in classifications:
                raise NlrSummaryError(
                    f"{args.article_matrix.name}:{line_number}: invalid unit/gene key"
                )
            classifications[key] = {
                "classification": row["manuscript_classification"],
                "refined_cause": row["refined_decayed_cause"],
                "evidence_level": row["refined_cause_evidence_level"],
            }
    if matrix_rows != args.expected_units * args.expected_reference_genes:
        raise NlrSummaryError(f"article matrix row count changed: {matrix_rows}")
    if len(classifications) != args.expected_units * args.expected_reference_nlrs:
        raise NlrSummaryError("article matrix does not cover every reference NLR and unit")

    unit_output: list[dict[str, object]] = []
    positive_output: list[dict[str, object]] = []
    type_output: list[dict[str, object]] = []
    universe_id = (
        "clem_scandens_article_nonshared_nlr_v1_"
        + hashlib.sha256(
            ("\n".join(sorted(nonshared_nlrs)) + "\n").encode("utf-8")
        ).hexdigest()[:12]
    )
    for metadata in metadata_rows:
        unit = metadata["assembly_unit_id"]
        counts: Counter[str] = Counter()
        type_counts: Counter[str] = Counter()
        for gene in sorted(nonshared_nlrs):
            evidence = classifications[(unit, gene)]
            classification = evidence["classification"]
            counts[classification] += 1
            if classification in POSITIVE:
                grouped_type = loss_type(evidence["refined_cause"])
                type_counts[grouped_type] += 1
                positive_output.append(
                    {
                        "assembly_unit_id": unit,
                        "reference_nlr_id": gene,
                        "primary_classification": classification,
                        "refined_cause": evidence["refined_cause"],
                        "refined_cause_evidence_level": evidence["evidence_level"],
                        "loss_type_group": grouped_type,
                        "reference_nlr_universe_id": universe_id,
                    }
                )
        unexpected = set(counts) - (RESOLVED | {"not_called_loss"})
        if unexpected:
            raise NlrSummaryError(f"{unit}: unexpected article classes {unexpected}")
        denominator = sum(counts[item] for item in RESOLVED)
        positive = sum(counts[item] for item in POSITIVE)
        if sum(type_counts.values()) != positive:
            raise NlrSummaryError(f"{unit}: loss-type counts do not close")
        percentage = 100.0 * positive / denominator if denominator else None
        unit_output.append(
            {
                "analysis_cohort": "article_method_23_units_nonshared_nlr_v1",
                "cohort_role": "primary",
                "assembly_unit_id": unit,
                "biological_species": metadata["biological_species"],
                "haplotype_or_subgenome": metadata["haplotype_or_subgenome"],
                "assembly_scope": metadata["assembly_scope"],
                "total_nlr_count": int(repertoire[unit]["total_nlr_count"]),
                "article_retained_reference_nlr_count": counts["retained"],
                "article_decayed_reference_nlr_loss_count": counts["decayed"],
                "article_deleted_reference_nlr_loss_count": counts["deleted"],
                "no_qualifying_translated_hit_count": type_counts[
                    "no_qualifying_translated_hit"
                ],
                "frameshift_supported_count": type_counts[
                    "frameshift_supported"
                ],
                "inframe_stop_supported_count": type_counts[
                    "inframe_stop_supported"
                ],
                "frameshift_and_stop_supported_count": type_counts[
                    "frameshift_and_stop_supported"
                ],
                "truncation_or_partial_alignment_candidate_count": type_counts[
                    "truncation_or_partial_alignment_candidate"
                ],
                "residual_sequence_mechanism_unresolved_count": type_counts[
                    "residual_sequence_mechanism_unresolved"
                ],
                "positive_reference_nlr_loss_count": positive,
                "callable_reference_nlr_denominator": denominator,
                "positive_reference_nlr_loss_percentage": (
                    f"{percentage:.6f}" if percentage is not None else ""
                ),
                "percentage_status": (
                    "defined" if percentage is not None else "undefined_zero_denominator"
                ),
                "reference_nlr_universe_id": universe_id,
            }
        )
        for grouped_type in LOSS_TYPE_ORDER:
            type_output.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": metadata["biological_species"],
                    "haplotype_or_subgenome": metadata[
                        "haplotype_or_subgenome"
                    ],
                    "loss_type_group": grouped_type,
                    "positive_reference_nlr_loss_count": type_counts[grouped_type],
                    "reference_nlr_universe_id": universe_id,
                }
            )

    universe_output = [
        {
            "reference_nlr_id": gene,
            "included_in_article_nonshared_analysis": (
                "false" if gene in shared_nlrs else "true"
            ),
            "exclusion_reason": (
                "positive_in_all_23_units_under_article_method"
                if gene in shared_nlrs
                else ""
            ),
            "reference_nlr_universe_id": universe_id,
        }
        for gene in sorted(reference_nlrs)
    ]

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        unit_path = staging / "article_nlr_unit_summary.tsv"
        universe_path = staging / "article_nlr_reference_universe.tsv"
        positive_path = staging / "positive_article_nlr_loss_calls.tsv.gz"
        type_path = staging / "nlr_loss_type_summary.tsv"
        write_tsv(unit_path, list(unit_output[0]), unit_output)
        write_tsv(universe_path, list(universe_output[0]), universe_output)
        write_tsv_gz(
            positive_path,
            list(positive_output[0]),
            positive_output,
        )
        write_tsv(type_path, list(type_output[0]), type_output)
        output_paths = [unit_path, universe_path, positive_path, type_path]
        all_type_counts = Counter(
            row["loss_type_group"] for row in positive_output
        )
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_ARTICLE_METHOD_NLR_SUMMARY",
            "loss_numerator": "article-method decayed + deleted",
            "loss_denominator": "retained + decayed + deleted",
            "not_called_policy": "excluded from numerator and denominator",
            "species_aggregation": "not performed",
            "assembly_units": args.expected_units,
            "biological_species": args.expected_species,
            "reference_nlrs": args.expected_reference_nlrs,
            "article_shared_reference_nlrs_excluded": len(shared_nlrs),
            "article_nonshared_reference_nlrs": len(nonshared_nlrs),
            "positive_unit_gene_calls": len(positive_output),
            "loss_type_groups": {
                grouped_type: all_type_counts[grouped_type]
                for grouped_type in LOSS_TYPE_ORDER
            },
            "callable_unit_gene_denominator": sum(
                int(row["callable_reference_nlr_denominator"])
                for row in unit_output
            ),
            "inputs": [
                {"role": role, "basename": path.name, "sha256": sha256(path)}
                for role, path in (
                    ("article_matrix", args.article_matrix),
                    ("article_shared_genes", args.shared_genes),
                    ("reference_nlr_universe", args.reference_nlr_universe),
                    ("repertoire_counts", args.repertoire_counts),
                    ("unit_metadata", args.unit_metadata),
                )
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in output_paths
            ],
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, args.output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (NlrSummaryError, OSError, csv.Error, UnicodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
