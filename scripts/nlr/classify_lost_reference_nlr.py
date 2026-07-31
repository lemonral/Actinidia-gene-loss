#!/usr/bin/env python3
"""Classify non-shared reference-NLR loss calls by NLR structural class."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path


NLR_CLASSES = (
    "CC-NBARC",
    "CC-NBARC-LRR",
    "NBARC",
    "NBARC-LRR",
    "TIR",
    "TIR-LRR",
    "TIR-NBARC",
    "TIR-NBARC-LRR",
    "TIR-CC-NBARC-LRR",
)
RESOLVED = {"retained", "decayed", "deleted"}
POSITIVE = {"decayed", "deleted"}
ALL_STATES = RESOLVED | {"not_called_loss"}


class NlrClassError(ValueError):
    """Raised when NLR structural classes and frozen loss calls do not close."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loss-matrix", required=True, type=Path)
    parser.add_argument("--reference-nlr-universe", required=True, type=Path)
    parser.add_argument("--reference-nlr-calls", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--nlr-output-root", required=True, type=Path)
    parser.add_argument("--nlr-unit-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-nlrs", type=int, default=214)
    parser.add_argument("--expected-nonshared-reference-nlrs", type=int, default=76)
    parser.add_argument("--expected-positive-calls", type=int, default=254)
    parser.add_argument("--expected-resolved-denominator", type=int, default=1738)
    parser.add_argument("--expected-repertoire-nlrs", type=int, default=6034)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise NlrClassError(f"missing or empty input: {path}")
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    binding(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = list(reader.fieldnames or [])
    if not rows or not fields or len(fields) != len(set(fields)):
        raise NlrClassError(f"invalid TSV: {path.name}")
    return rows, fields


def require(fields: list[str], required: set[str], label: str) -> None:
    missing = sorted(required - set(fields))
    if missing:
        raise NlrClassError(f"{label} missing fields: {', '.join(missing)}")


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


def read_reference_classes(path: Path) -> dict[str, str]:
    binding(path)
    classes: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise NlrClassError(
                    f"{path.name}:{line_number}: expected seven fields"
                )
            gene, nlr_class = fields[0], fields[2]
            if nlr_class not in NLR_CLASSES:
                raise NlrClassError(
                    f"{path.name}:{line_number}: unknown NLR class {nlr_class!r}"
                )
            if gene in classes:
                raise NlrClassError(
                    f"{path.name}:{line_number}: duplicate reference NLR {gene}"
                )
            classes[gene] = nlr_class
    if not classes:
        raise NlrClassError("reference NLR calls are empty")
    return classes


def read_nlr_class_counts(path: Path) -> Counter[str]:
    binding(path)
    counts: Counter[str] = Counter()
    nlr_ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise NlrClassError(
                    f"{path.name}:{line_number}: expected seven fields"
                )
            nlr_id, nlr_class = fields[1], fields[2]
            if nlr_class not in NLR_CLASSES:
                raise NlrClassError(
                    f"{path.name}:{line_number}: unknown NLR class {nlr_class!r}"
                )
            if nlr_id in nlr_ids:
                raise NlrClassError(
                    f"{path.name}:{line_number}: duplicate NLR call {nlr_id}"
                )
            nlr_ids.add(nlr_id)
            counts[nlr_class] += 1
    if not nlr_ids:
        raise NlrClassError(f"empty NLR repertoire: {path}")
    return counts


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise NlrClassError(
            f"refusing to overwrite output directory: {args.output_dir}"
        )
    input_paths = [
        args.loss_matrix.resolve(),
        args.reference_nlr_universe.resolve(),
        args.reference_nlr_calls.resolve(),
        args.unit_metadata.resolve(),
        args.nlr_unit_summary.resolve(),
    ]
    for path in input_paths:
        binding(path)
    nlr_output_root = args.nlr_output_root.resolve()
    if not nlr_output_root.is_dir():
        raise NlrClassError(
            f"missing NLR-Annotator output root: {nlr_output_root}"
        )

    metadata_rows, metadata_fields = read_tsv(args.unit_metadata)
    require(
        metadata_fields,
        {
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "include",
        },
        args.unit_metadata.name,
    )
    metadata_rows = [
        row for row in metadata_rows if row["include"].lower() == "true"
    ]
    units = [row["assembly_unit_id"] for row in metadata_rows]
    metadata = {row["assembly_unit_id"]: row for row in metadata_rows}
    if len(units) != args.expected_units or len(metadata) != len(units):
        raise NlrClassError("unit metadata does not contain the exact cohort")

    universe_rows, universe_fields = read_tsv(args.reference_nlr_universe)
    require(
        universe_fields,
        {
            "reference_nlr_id",
            "included_in_article_nonshared_analysis",
            "exclusion_reason",
        },
        args.reference_nlr_universe.name,
    )
    universe = {row["reference_nlr_id"]: row for row in universe_rows}
    if len(universe) != args.expected_reference_nlrs:
        raise NlrClassError("reference NLR universe count changed")
    nonshared = {
        gene
        for gene, row in universe.items()
        if row["included_in_article_nonshared_analysis"].lower() == "true"
    }
    if len(nonshared) != args.expected_nonshared_reference_nlrs:
        raise NlrClassError("non-shared reference NLR count changed")

    reference_classes = read_reference_classes(args.reference_nlr_calls)
    if set(reference_classes) != set(universe):
        raise NlrClassError(
            "reference NLR structural classes do not match the exact universe"
        )

    unit_summary_rows, unit_summary_fields = read_tsv(args.nlr_unit_summary)
    require(
        unit_summary_fields,
        {"assembly_unit_id", "total_nlr_count"},
        args.nlr_unit_summary.name,
    )
    declared_repertoire_totals = {
        row["assembly_unit_id"]: int(row["total_nlr_count"])
        for row in unit_summary_rows
    }
    if (
        set(declared_repertoire_totals) != set(units)
        or len(declared_repertoire_totals) != len(unit_summary_rows)
    ):
        raise NlrClassError("NLR unit summary does not match the exact cohort")

    repertoire_counts: dict[str, Counter[str]] = {}
    repertoire_bindings: list[dict[str, object]] = []
    for unit in units:
        calls_path = nlr_output_root / unit / "nlr_calls.txt"
        counts = read_nlr_class_counts(calls_path)
        observed_total = sum(counts.values())
        if observed_total != declared_repertoire_totals[unit]:
            raise NlrClassError(
                f"{unit}: NLR repertoire total {observed_total} does not match "
                f"declared total {declared_repertoire_totals[unit]}"
            )
        repertoire_counts[unit] = counts
        repertoire_bindings.append(
            {
                "assembly_unit_id": unit,
                **binding(calls_path),
            }
        )
    repertoire_total = sum(
        sum(counts.values()) for counts in repertoire_counts.values()
    )
    if repertoire_total != args.expected_repertoire_nlrs:
        raise NlrClassError(
            f"total NLR repertoire count changed: {repertoire_total}"
        )

    classifications: dict[tuple[str, str], dict[str, str]] = {}
    with open_text(args.loss_matrix) as handle:
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
            args.loss_matrix.name,
        )
        for line_number, row in enumerate(reader, 2):
            gene = row["reference_gene_id"]
            if gene not in nonshared:
                continue
            unit = row["assembly_unit_id"]
            key = (unit, gene)
            classification = row["manuscript_classification"]
            if (
                unit not in metadata
                or key in classifications
                or classification not in ALL_STATES
            ):
                raise NlrClassError(
                    f"{args.loss_matrix.name}:{line_number}: invalid loss row"
                )
            classifications[key] = {
                "classification": classification,
                "refined_cause": row["refined_decayed_cause"],
                "refined_cause_evidence_level": row[
                    "refined_cause_evidence_level"
                ],
            }
    expected_pairs = args.expected_units * args.expected_nonshared_reference_nlrs
    if len(classifications) != expected_pairs:
        raise NlrClassError(
            f"loss matrix covers {len(classifications)} rather than "
            f"{expected_pairs} unit-gene pairs"
        )

    class_nonshared_counts = Counter(
        reference_classes[gene] for gene in nonshared
    )
    class_all_counts = Counter(reference_classes.values())
    summary_rows: list[dict[str, object]] = []
    positive_rows: list[dict[str, object]] = []
    resolved_total = 0
    for unit in units:
        for nlr_class in NLR_CLASSES:
            genes = sorted(
                gene
                for gene in nonshared
                if reference_classes[gene] == nlr_class
            )
            counts: Counter[str] = Counter()
            for gene in genes:
                evidence = classifications[(unit, gene)]
                classification = evidence["classification"]
                counts[classification] += 1
                if classification in POSITIVE:
                    positive_rows.append(
                        {
                            "assembly_unit_id": unit,
                            "biological_species": metadata[unit][
                                "biological_species"
                            ],
                            "haplotype_or_subgenome": metadata[unit][
                                "haplotype_or_subgenome"
                            ],
                            "reference_nlr_id": gene,
                            "reference_nlr_class": nlr_class,
                            "primary_classification": classification,
                            "refined_cause": evidence["refined_cause"],
                            "refined_cause_evidence_level": evidence[
                                "refined_cause_evidence_level"
                            ],
                        }
                    )
            resolved = sum(counts[state] for state in RESOLVED)
            positive = sum(counts[state] for state in POSITIVE)
            resolved_total += resolved
            percentage = 100.0 * positive / resolved if resolved else None
            summary_rows.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": metadata[unit]["biological_species"],
                    "haplotype_or_subgenome": metadata[unit][
                        "haplotype_or_subgenome"
                    ],
                    "reference_nlr_class": nlr_class,
                    "all_reference_nlr_genes_in_class": class_all_counts[nlr_class],
                    "nonshared_reference_nlr_genes_in_class": class_nonshared_counts[
                        nlr_class
                    ],
                    "resolved_unit_gene_denominator": resolved,
                    "retained_count": counts["retained"],
                    "decayed_loss_count": counts["decayed"],
                    "deleted_loss_count": counts["deleted"],
                    "positive_loss_count": positive,
                    "not_called_count": counts["not_called_loss"],
                    "positive_loss_percentage": (
                        f"{percentage:.6f}" if percentage is not None else ""
                    ),
                }
            )
    if len(positive_rows) != args.expected_positive_calls:
        raise NlrClassError(
            f"positive NLR loss count changed: {len(positive_rows)}"
        )
    if resolved_total != args.expected_resolved_denominator:
        raise NlrClassError(
            f"resolved NLR denominator changed: {resolved_total}"
        )

    class_rows = [
        {
            "reference_nlr_id": gene,
            "reference_nlr_class": reference_classes[gene],
            "included_in_nonshared_loss_analysis": str(gene in nonshared).lower(),
            "exclusion_reason": universe[gene]["exclusion_reason"],
        }
        for gene in sorted(universe)
    ]
    shared_genes = set(universe) - nonshared
    if any(
        universe[gene]["exclusion_reason"]
        != "positive_in_all_23_units_under_article_method"
        for gene in shared_genes
    ):
        raise NlrClassError("shared NLR exclusion reason changed")
    shared_counts = Counter(reference_classes[gene] for gene in shared_genes)
    shared_rows = [
        {
            "cohort": "shared_positive_all_23_units",
            "display_label": "Shared loss",
            "reference_nlr_class": nlr_class,
            "shared_reference_nlr_gene_count": shared_counts[nlr_class],
        }
        for nlr_class in NLR_CLASSES
    ]
    repertoire_rows = [
        {
            "assembly_unit_id": unit,
            "biological_species": metadata[unit]["biological_species"],
            "haplotype_or_subgenome": metadata[unit][
                "haplotype_or_subgenome"
            ],
            "reference_nlr_class": nlr_class,
            "nlr_gene_count": repertoire_counts[unit][nlr_class],
            "total_nlr_count": declared_repertoire_totals[unit],
        }
        for unit in units
        for nlr_class in NLR_CLASSES
    ]
    summary_index = {
        (row["assembly_unit_id"], row["reference_nlr_class"]): row
        for row in summary_rows
    }
    rate_rows: list[dict[str, object]] = []
    for unit in units:
        for nlr_class in NLR_CLASSES:
            row = summary_index[(unit, nlr_class)]
            shared_count = shared_counts[nlr_class]
            nonshared_positive = int(row["positive_loss_count"])
            nonshared_resolved = int(row["resolved_unit_gene_denominator"])
            all_positive = shared_count + nonshared_positive
            all_resolved = shared_count + nonshared_resolved
            if all_positive > all_resolved:
                raise NlrClassError(
                    f"{unit}/{nlr_class}: positive NLR calls exceed "
                    "resolved opportunities"
                )
            all_percentage = (
                100.0 * all_positive / all_resolved
                if all_resolved
                else None
            )
            nonshared_percentage = (
                100.0 * nonshared_positive / nonshared_resolved
                if nonshared_resolved
                else None
            )
            rate_rows.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": metadata[unit][
                        "biological_species"
                    ],
                    "haplotype_or_subgenome": metadata[unit][
                        "haplotype_or_subgenome"
                    ],
                    "reference_nlr_class": nlr_class,
                    "shared_loss_count": shared_count,
                    "nonshared_positive_loss_count": nonshared_positive,
                    "nonshared_resolved_denominator": nonshared_resolved,
                    "all_positive_loss_count": all_positive,
                    "all_resolved_denominator": all_resolved,
                    "all_loss_percentage": (
                        f"{all_percentage:.6f}"
                        if all_percentage is not None
                        else ""
                    ),
                    "nonshared_loss_percentage": (
                        f"{nonshared_percentage:.6f}"
                        if nonshared_percentage is not None
                        else ""
                    ),
                }
            )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        summary_path = staging / "lost_nlr_structural_class_summary.tsv"
        positive_path = staging / "lost_nlr_structural_class_calls.tsv"
        classes_path = staging / "reference_nlr_structural_classes.tsv"
        shared_path = staging / "shared_nlr_structural_class_summary.tsv"
        repertoire_path = (
            staging / "nlr_repertoire_structural_class_summary.tsv"
        )
        rate_path = staging / "nlr_structural_class_loss_rates.tsv"
        write_tsv(summary_path, list(summary_rows[0]), summary_rows)
        write_tsv(positive_path, list(positive_rows[0]), positive_rows)
        write_tsv(classes_path, list(class_rows[0]), class_rows)
        write_tsv(shared_path, list(shared_rows[0]), shared_rows)
        write_tsv(
            repertoire_path,
            list(repertoire_rows[0]),
            repertoire_rows,
        )
        write_tsv(rate_path, list(rate_rows[0]), rate_rows)
        output_paths = [
            summary_path,
            positive_path,
            classes_path,
            shared_path,
            repertoire_path,
            rate_path,
        ]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_LOST_NLR_STRUCTURAL_CLASSIFICATION",
            "assembly_units": len(units),
            "reference_nlrs": len(universe),
            "nonshared_reference_nlrs": len(nonshared),
            "shared_reference_nlrs_excluded": len(universe) - len(nonshared),
            "shared_reference_nlr_class_counts": {
                nlr_class: shared_counts[nlr_class]
                for nlr_class in NLR_CLASSES
            },
            "positive_unit_gene_calls": len(positive_rows),
            "resolved_unit_gene_denominator": resolved_total,
            "repertoire_nlr_calls": repertoire_total,
            "nlr_structural_classes": list(NLR_CLASSES),
            "reference_nlr_class_counts": {
                nlr_class: class_all_counts[nlr_class]
                for nlr_class in NLR_CLASSES
            },
            "nonshared_reference_nlr_class_counts": {
                nlr_class: class_nonshared_counts[nlr_class]
                for nlr_class in NLR_CLASSES
            },
            "loss_numerator": "decayed + deleted",
            "loss_denominator": "retained + decayed + deleted within NLR class",
            "all_loss_rate_denominator": (
                "shared positive reference NLRs + resolved non-shared "
                "reference NLR opportunities within class"
            ),
            "nonshared_loss_rate_denominator": (
                "resolved non-shared reference NLR opportunities within class"
            ),
            "not_called_policy": "excluded from numerator and denominator",
            "species_aggregation": "not performed",
            "inputs": [binding(path) for path in input_paths],
            "repertoire_inputs": repertoire_bindings,
            "outputs": [binding(path) for path in output_paths],
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
    except (NlrClassError, OSError, UnicodeError, csv.Error) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
