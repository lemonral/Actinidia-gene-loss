#!/usr/bin/env python3
"""Summarize article-method losses and orthogonal strict disruption evidence.

Article-method ``decayed + deleted`` remains the publication trend numerator.
Strict pseudogenized calls are an evidence refinement and are never added to
that numerator a second time.  Branch rows use exact article-method topology
events and report strict disruption support separately.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


csv.field_size_limit(sys.maxsize)


class EvidenceError(ValueError):
    """Raised when frozen article and strict-evidence tables disagree."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article-unit-summary", required=True, type=Path)
    parser.add_argument("--uniform-unit-summary", required=True, type=Path)
    parser.add_argument("--classification-crosswalk", required=True, type=Path)
    parser.add_argument("--refined-cause-summary", required=True, type=Path)
    parser.add_argument("--strict-disruption-calls", required=True, type=Path)
    parser.add_argument("--tree-loss-events", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-strict-calls", type=int, default=20046)
    parser.add_argument("--expected-branch-events", type=int, default=26729)
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
        raise EvidenceError(f"missing or empty input: {path.name}")
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    if not fields or not rows or len(fields) != len(set(fields)):
        raise EvidenceError(f"invalid TSV input: {path.name}")
    return rows, fields


def require(fields: list[str], needed: set[str], label: str) -> None:
    missing = sorted(needed - set(fields))
    if missing:
        raise EvidenceError(f"{label} missing columns: {', '.join(missing)}")


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


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise EvidenceError(f"output directory already exists: {args.output_dir}")

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
        row for row in metadata_rows if row["include"].strip().lower() == "true"
    ]
    metadata = {row["assembly_unit_id"]: row for row in metadata_rows}
    if len(metadata_rows) != args.expected_units or len(metadata) != args.expected_units:
        raise EvidenceError("unit metadata is not the exact unique 23-unit cohort")
    species_units: dict[str, set[str]] = defaultdict(set)
    for row in metadata_rows:
        species_units[row["biological_species"]].add(row["assembly_unit_id"])
    if len(species_units) != 13:
        raise EvidenceError("unit metadata does not close to 13 lineages")

    article_rows, article_fields = read_tsv(args.article_unit_summary)
    require(
        article_fields,
        {
            "assembly_unit_id",
            "retained",
            "decayed",
            "deleted",
            "not_called_loss",
            "manuscript_positive_loss",
            "manuscript_resolved_denominator",
        },
        args.article_unit_summary.name,
    )
    article = {row["assembly_unit_id"]: row for row in article_rows}
    if len(article_rows) != args.expected_units or set(article) != set(metadata):
        raise EvidenceError("article summary unit universe differs from metadata")

    uniform_rows, uniform_fields = read_tsv(args.uniform_unit_summary)
    require(
        uniform_fields,
        {
            "assembly_unit_id",
            "retained",
            "deleted",
            "pseudogenized",
            "uncertain",
        },
        args.uniform_unit_summary.name,
    )
    uniform = {row["assembly_unit_id"]: row for row in uniform_rows}
    if len(uniform_rows) != args.expected_units or set(uniform) != set(metadata):
        raise EvidenceError("uniform summary unit universe differs from metadata")

    strict_rows, strict_fields = read_tsv(args.strict_disruption_calls)
    require(
        strict_fields,
        {
            "reference_gene_id",
            "assembly_unit_id",
            "frameshift_events",
            "inframe_stop_codons",
        },
        args.strict_disruption_calls.name,
    )
    if len(strict_rows) != args.expected_strict_calls:
        raise EvidenceError(
            f"strict disruption rows={len(strict_rows)}; expected {args.expected_strict_calls}"
        )
    strict_by_gene: dict[str, set[str]] = defaultdict(set)
    strict_by_unit: Counter[str] = Counter()
    strict_types: Counter[str] = Counter()
    seen_strict: set[tuple[str, str]] = set()
    for line_number, row in enumerate(strict_rows, 2):
        gene = row["reference_gene_id"].strip()
        unit = row["assembly_unit_id"].strip()
        key = (gene, unit)
        if not gene or unit not in metadata or key in seen_strict:
            raise EvidenceError(
                f"{args.strict_disruption_calls.name}:{line_number}: invalid strict call"
            )
        seen_strict.add(key)
        try:
            frameshifts = int(row["frameshift_events"])
            stops = int(row["inframe_stop_codons"])
        except ValueError as error:
            raise EvidenceError(
                f"{args.strict_disruption_calls.name}:{line_number}: invalid event counts"
            ) from error
        if frameshifts < 0 or stops < 0 or frameshifts + stops < 1:
            raise EvidenceError(
                f"{args.strict_disruption_calls.name}:{line_number}: unsupported strict event"
            )
        if frameshifts and stops:
            evidence_type = "frameshift_and_inframe_stop"
        elif frameshifts:
            evidence_type = "frameshift_only"
        else:
            evidence_type = "inframe_stop_only"
        strict_types[evidence_type] += 1
        strict_by_gene[gene].add(unit)
        strict_by_unit[unit] += 1
    if sum(strict_types.values()) != args.expected_strict_calls:
        raise EvidenceError("strict evidence types do not close")

    unit_output: list[dict[str, object]] = []
    article_totals: Counter[str] = Counter()
    for unit in metadata:
        row = article[unit]
        values = {
            key: int(row[key])
            for key in (
                "retained",
                "decayed",
                "deleted",
                "not_called_loss",
                "manuscript_positive_loss",
                "manuscript_resolved_denominator",
            )
        }
        if (
            values["retained"]
            + values["decayed"]
            + values["deleted"]
            + values["not_called_loss"]
            != args.expected_reference_genes
            or values["decayed"] + values["deleted"]
            != values["manuscript_positive_loss"]
            or values["retained"] + values["manuscript_positive_loss"]
            != values["manuscript_resolved_denominator"]
        ):
            raise EvidenceError(f"{unit}: article counts do not close")
        strict_count = strict_by_unit[unit]
        if strict_count != int(uniform[unit]["pseudogenized"]):
            raise EvidenceError(f"{unit}: strict call count differs from uniform summary")
        article_totals.update(
            {
                "retained": values["retained"],
                "decayed": values["decayed"],
                "deleted": values["deleted"],
                "not_called_loss": values["not_called_loss"],
            }
        )
        unit_output.append(
            {
                "assembly_unit_id": unit,
                "biological_species": metadata[unit]["biological_species"],
                "haplotype_or_subgenome": metadata[unit][
                    "haplotype_or_subgenome"
                ],
                "retained": values["retained"],
                "article_decayed": values["decayed"],
                "article_deleted": values["deleted"],
                "article_positive_loss": values["manuscript_positive_loss"],
                "article_resolved_denominator": values[
                    "manuscript_resolved_denominator"
                ],
                "article_loss_rate": (
                    values["manuscript_positive_loss"]
                    / values["manuscript_resolved_denominator"]
                ),
                "strict_pseudogenized_evidence": strict_count,
                "strict_evidence_fraction_of_article_positive": (
                    strict_count / values["manuscript_positive_loss"]
                ),
                "not_called_loss": values["not_called_loss"],
            }
        )
    if article_totals != Counter(
        {
            "retained": 633957,
            "decayed": 171866,
            "deleted": 7961,
            "not_called_loss": 3797,
        }
    ):
        raise EvidenceError(f"article totals changed: {dict(article_totals)}")

    crosswalk_rows, crosswalk_fields = read_tsv(args.classification_crosswalk)
    require(
        crosswalk_fields,
        {
            "manuscript_classification",
            "uniform_classification",
            "refined_decayed_cause",
            "row_count",
        },
        args.classification_crosswalk.name,
    )
    strict_article_positive = 0
    strict_outside_article_positive = 0
    for row in crosswalk_rows:
        count = int(row["row_count"])
        if row["uniform_classification"] != "pseudogenized":
            continue
        if row["manuscript_classification"] in {"decayed", "deleted"}:
            strict_article_positive += count
        else:
            strict_outside_article_positive += count
    if (
        strict_article_positive + strict_outside_article_positive
        != args.expected_strict_calls
        or strict_outside_article_positive != 1
    ):
        raise EvidenceError("strict/article crosswalk does not close")

    refined_rows, refined_fields = read_tsv(args.refined_cause_summary)
    require(
        refined_fields,
        {"assembly_unit_id", "refined_decayed_cause", "row_count"},
        args.refined_cause_summary.name,
    )
    refined: Counter[str] = Counter()
    refined_units: set[str] = set()
    for row in refined_rows:
        if row["assembly_unit_id"] not in metadata:
            raise EvidenceError("refined cause table contains an unknown unit")
        refined_units.add(row["assembly_unit_id"])
        refined[row["refined_decayed_cause"]] += int(row["row_count"])
    if refined_units != set(metadata) or sum(refined.values()) != (
        args.expected_units * args.expected_reference_genes
    ):
        raise EvidenceError("refined cause table does not close to the complete grid")
    refined_output = [
        {"refined_cause": key, "unit_gene_rows": value}
        for key, value in sorted(refined.items(), key=lambda item: (-item[1], item[0]))
    ]

    event_rows, event_fields = read_tsv(args.tree_loss_events)
    require(
        event_fields,
        {
            "reference_gene_id",
            "branch_id",
            "branch_type",
            "descendant_lineage_count",
            "descendant_lineages",
        },
        args.tree_loss_events.name,
    )
    if len(event_rows) != args.expected_branch_events:
        raise EvidenceError(
            f"tree event rows={len(event_rows)}; expected {args.expected_branch_events}"
        )
    branch_genes: dict[str, set[str]] = defaultdict(set)
    branch_meta: dict[str, tuple[str, tuple[str, ...]]] = {}
    for line_number, row in enumerate(event_rows, 2):
        gene = row["reference_gene_id"].strip()
        branch = row["branch_id"].strip()
        branch_type = row["branch_type"].strip()
        descendants = tuple(row["descendant_lineages"].split(";"))
        if (
            not gene
            or not branch
            or branch_type not in {"terminal", "internal"}
            or len(descendants) != int(row["descendant_lineage_count"])
            or any(species not in species_units for species in descendants)
        ):
            raise EvidenceError(
                f"{args.tree_loss_events.name}:{line_number}: invalid branch event"
            )
        if gene in branch_genes[branch]:
            raise EvidenceError(
                f"{args.tree_loss_events.name}:{line_number}: duplicate gene-branch event"
            )
        branch_genes[branch].add(gene)
        current_meta = (branch_type, descendants)
        if branch in branch_meta and branch_meta[branch] != current_meta:
            raise EvidenceError("branch metadata changed between event rows")
        branch_meta[branch] = current_meta
    if sum(map(len, branch_genes.values())) != args.expected_branch_events:
        raise EvidenceError("unique tree events do not close")

    branch_output: list[dict[str, object]] = []
    for branch in sorted(
        branch_genes,
        key=lambda value: (
            0 if branch_meta[value][0] == "terminal" else 1,
            branch_meta[value][1],
        ),
    ):
        branch_type, descendants = branch_meta[branch]
        any_supported = 0
        all_supported = 0
        for gene in branch_genes[branch]:
            supporting_units = strict_by_gene.get(gene, set())
            lineage_support = [
                bool(species_units[species].intersection(supporting_units))
                for species in descendants
            ]
            any_supported += int(any(lineage_support))
            all_supported += int(all(lineage_support))
        article_gene_count = len(branch_genes[branch])
        branch_output.append(
            {
                "branch_id": branch,
                "branch_type": branch_type,
                "descendant_lineage_count": len(descendants),
                "descendant_lineages": ";".join(descendants),
                "article_method_loss_gene_count": article_gene_count,
                "strict_supported_any_descendant_lineage": any_supported,
                "strict_supported_all_descendant_lineages": all_supported,
                "strict_all_lineages_fraction": all_supported
                / article_gene_count,
                "strict_support_definition": (
                    "at_least_one_strict_pseudogenized_unit_per_descendant_lineage"
                ),
            }
        )
    if Counter(row["branch_type"] for row in branch_output) != Counter(
        {"terminal": 13, "internal": 12}
    ):
        raise EvidenceError("branch output is not 13 terminal plus 12 internal rows")

    strict_output = [
        {
            "strict_disruption_type": evidence_type,
            "strict_pseudogenized_unit_gene_rows": strict_types[evidence_type],
        }
        for evidence_type in (
            "frameshift_only",
            "inframe_stop_only",
            "frameshift_and_inframe_stop",
        )
    ]

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        unit_path = staging / "unit_loss_evidence_summary.tsv"
        branch_path = staging / "branch_loss_evidence_summary.tsv"
        strict_path = staging / "strict_disruption_type_summary.tsv"
        refined_path = staging / "refined_cause_summary.tsv"
        write_tsv(unit_path, list(unit_output[0]), unit_output)
        write_tsv(branch_path, list(branch_output[0]), branch_output)
        write_tsv(strict_path, list(strict_output[0]), strict_output)
        write_tsv(refined_path, list(refined_output[0]), refined_output)
        output_paths = [unit_path, branch_path, strict_path, refined_path]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_ARTICLE_LOSS_WITH_STRICT_EVIDENCE_SUMMARY",
            "article_loss_numerator": "decayed + deleted",
            "strict_evidence_role": (
                "orthogonal high-confidence disruption support; never additive "
                "to the article numerator"
            ),
            "article_unit_gene_rows": args.expected_units
            * args.expected_reference_genes,
            "article_positive_unit_gene_rows": article_totals["decayed"]
            + article_totals["deleted"],
            "strict_pseudogenized_unit_gene_rows": args.expected_strict_calls,
            "strict_rows_within_article_positive": strict_article_positive,
            "strict_rows_outside_article_positive": strict_outside_article_positive,
            "tree_branch_event_rows": args.expected_branch_events,
            "terminal_branches": 13,
            "internal_branches": 12,
            "inputs": [
                {"role": role, "basename": path.name, "sha256": sha256(path)}
                for role, path in (
                    ("article_unit_summary", args.article_unit_summary),
                    ("uniform_unit_summary", args.uniform_unit_summary),
                    ("classification_crosswalk", args.classification_crosswalk),
                    ("refined_cause_summary", args.refined_cause_summary),
                    ("strict_disruption_calls", args.strict_disruption_calls),
                    ("tree_loss_events", args.tree_loss_events),
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
    except (EvidenceError, OSError, csv.Error, UnicodeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
