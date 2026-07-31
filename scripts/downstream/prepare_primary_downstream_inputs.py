#!/usr/bin/env python3
"""Prepare resolved-loss, expression, and copy-number inputs atomically.

The primary loss matrix is a complete reference-gene by assembly-unit grid,
but ``uncertain`` calls are not valid rate-denominator observations.  This
adapter therefore publishes only resolved ``deleted``, ``pseudogenized``, and ``retained`` rows
for expression/copy-number modelling while retaining per-unit uncertainty
counts in an audit table.  It also extracts one explicitly named featureCounts
sample as raw counts and freezes reference proteins absent from the historical
CD-HIT 0.90 membership file as an exclusion ledger.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd


LOSS_FIELDS = (
    "reference_gene_id",
    "assembly_unit_id",
    "classification",
    "callable",
    "evidence_source",
    "primary_search_state",
)
LEDGER_FIELDS = (
    "sample_id",
    "species",
    "ploidy",
    "analysis_role",
    "input_scope",
    "source_fasta",
    "output_basename",
    "expected_fasta_records",
)
PLOIDY_LABELS = {"2x": "diploid", "4x": "tetraploid", "6x": "hexaploid"}
CLUSTER_MEMBER = re.compile(r">(.+?)\.\.\.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loss-matrix", type=Path, required=True)
    parser.add_argument("--unit-ledger", type=Path, required=True)
    parser.add_argument("--reference-cds", type=Path, required=True)
    parser.add_argument("--shared-positive-genes", type=Path, required=True)
    parser.add_argument("--feature-counts", type=Path, required=True)
    parser.add_argument("--expression-sample-column", required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()
    current_has_bases = False
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(b">"):
                if identifiers and not current_has_bases:
                    raise ValueError(f"Empty FASTA record before line {line_number}")
                header = raw[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                try:
                    identifier = header.split()[0].decode("utf-8")
                except UnicodeDecodeError as exc:
                    raise ValueError(f"Non-UTF-8 FASTA header at line {line_number}") from exc
                if identifier in seen:
                    raise ValueError(f"Duplicate reference FASTA ID: {identifier}")
                seen.add(identifier)
                identifiers.append(identifier)
                current_has_bases = False
            elif raw.strip():
                if not identifiers:
                    raise ValueError("Sequence data precedes the first FASTA header")
                current_has_bases = True
    if identifiers and not current_has_bases:
        raise ValueError("Final FASTA record is empty")
    if not identifiers:
        raise ValueError("Reference CDS FASTA is empty")
    return identifiers


def read_unit_ledger(path: Path, expected_units: int) -> pd.DataFrame:
    ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if tuple(ledger.columns) != LEDGER_FIELDS:
        raise ValueError(f"Unexpected unit-ledger columns: {list(ledger.columns)}")
    targets = ledger.loc[ledger["analysis_role"] == "target_repertoire"].copy()
    if len(targets) != expected_units:
        raise ValueError(f"Observed {len(targets)} target units; expected {expected_units}")
    if targets["sample_id"].duplicated().any():
        raise ValueError("Unit ledger contains duplicate target sample_id values")
    if set(targets["input_scope"]) != {"whole_genome"}:
        raise ValueError("Every target unit must use whole_genome scope")
    unknown = sorted(set(targets["ploidy"]) - set(PLOIDY_LABELS))
    if unknown:
        raise ValueError(f"Unsupported ploidy labels: {unknown}")
    targets["ploidy_class"] = targets["ploidy"].map(PLOIDY_LABELS)
    return targets[["sample_id", "species", "ploidy", "ploidy_class"]]


def read_cluster_members(path: Path) -> set[str]:
    members: set[str] = set()
    current_cluster: str | None = None
    cluster_count = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current_cluster = line[1:]
                cluster_count += 1
                continue
            match = CLUSTER_MEMBER.search(line)
            if current_cluster is None or match is None:
                raise ValueError(f"Malformed CD-HIT record at line {line_number}")
            gene = match.group(1)
            if gene in members:
                raise ValueError(f"Reference protein appears in multiple clusters: {gene}")
            members.add(gene)
    if not members or cluster_count == 0:
        raise ValueError("No CD-HIT clusters were parsed")
    return members


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def main() -> int:
    args = parse_args()
    if args.expected_units < 1 or args.expected_reference_genes < 1:
        print("ERROR: expected counts must be positive", file=sys.stderr)
        return 2
    if args.output_dir.exists():
        print(f"ERROR: output directory already exists: {args.output_dir}", file=sys.stderr)
        return 2
    inputs = [
        args.loss_matrix,
        args.unit_ledger,
        args.reference_cds,
        args.shared_positive_genes,
        args.feature_counts,
        args.clusters,
    ]
    try:
        missing = [str(path) for path in inputs if not path.is_file()]
        if missing:
            raise ValueError(f"Missing input files: {missing}")
        reference_order = fasta_ids(args.reference_cds)
        if len(reference_order) != args.expected_reference_genes:
            raise ValueError(
                f"Observed {len(reference_order)} reference genes; "
                f"expected {args.expected_reference_genes}"
            )
        reference = set(reference_order)
        shared_table = pd.read_csv(
            args.shared_positive_genes, sep="\t", dtype=str, keep_default_na=False
        )
        if "reference_gene_id" not in shared_table.columns:
            raise ValueError("Shared-positive table lacks reference_gene_id")
        if shared_table["reference_gene_id"].duplicated().any():
            raise ValueError("Shared-positive table contains duplicate reference_gene_id values")
        shared_positive = set(shared_table["reference_gene_id"])
        if not shared_positive or not shared_positive.issubset(reference):
            raise ValueError("Shared-positive gene universe must be a nonempty reference subset")
        analysis_reference = reference - shared_positive
        analysis_order = [gene for gene in reference_order if gene in analysis_reference]
        units = read_unit_ledger(args.unit_ledger, args.expected_units)
        unit_ids = set(units["sample_id"])

        losses = pd.read_csv(
            args.loss_matrix, sep="\t", dtype=str, keep_default_na=False
        )
        missing_loss_fields = set(LOSS_FIELDS[:5]) - set(losses.columns)
        if missing_loss_fields:
            raise ValueError(f"Loss matrix lacks required columns: {sorted(missing_loss_fields)}")
        if "primary_search_state" not in losses.columns:
            if "evidence_reason" not in losses.columns:
                raise ValueError("Loss matrix lacks primary_search_state/evidence_reason")
            losses["primary_search_state"] = losses["evidence_reason"]
        expected_rows = args.expected_units * args.expected_reference_genes
        if len(losses) != expected_rows:
            raise ValueError(f"Observed {len(losses)} loss rows; expected {expected_rows}")
        if losses.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
            raise ValueError("Loss matrix contains duplicate reference-gene by unit rows")
        if set(losses["reference_gene_id"]) != reference:
            raise ValueError("Loss-matrix reference-gene universe does not match the CDS FASTA")
        if set(losses["assembly_unit_id"]) != unit_ids:
            raise ValueError("Loss-matrix assembly-unit universe does not match the unit ledger")
        allowed_classes = {"deleted", "pseudogenized", "retained", "uncertain"}
        if not set(losses["classification"]).issubset(allowed_classes):
            raise ValueError(
                f"Unexpected primary classification set: {sorted(set(losses['classification']))}"
            )
        nonshared_losses = losses.loc[
            ~losses["reference_gene_id"].isin(shared_positive)
        ].copy()
        resolved_mask = nonshared_losses["classification"].isin(
            {"deleted", "pseudogenized", "retained"}
        )
        if (nonshared_losses.loc[resolved_mask, "callable"] != "true").any():
            raise ValueError("Every resolved primary loss call must be callable")
        resolved = nonshared_losses.loc[resolved_mask].merge(
            units, left_on="assembly_unit_id", right_on="sample_id", validate="many_to_one"
        )
        resolved = resolved[
            [
                "reference_gene_id",
                "assembly_unit_id",
                "species",
                "ploidy_class",
                "ploidy",
                "classification",
                "callable",
                "evidence_source",
                "primary_search_state",
            ]
        ].rename(columns={"ploidy_class": "ploidy", "ploidy": "ploidy_level"})
        resolved = resolved.sort_values(["assembly_unit_id", "reference_gene_id"])

        audit = (
            nonshared_losses.groupby(["assembly_unit_id", "classification"]).size().unstack(fill_value=0)
            .reindex(columns=["deleted", "pseudogenized", "retained", "uncertain"], fill_value=0)
            .reset_index()
            .merge(units, left_on="assembly_unit_id", right_on="sample_id", validate="one_to_one")
        )
        audit["positive_loss"] = audit["deleted"] + audit["pseudogenized"]
        audit["resolved_denominator"] = audit["positive_loss"] + audit["retained"]
        audit["positive_loss_rate"] = audit["positive_loss"] / audit["resolved_denominator"]
        audit = audit[
            [
                "assembly_unit_id", "species", "ploidy", "ploidy_class",
                "deleted", "pseudogenized", "retained", "uncertain", "positive_loss", "resolved_denominator",
                "positive_loss_rate",
            ]
        ].sort_values("assembly_unit_id")

        counts = pd.read_csv(
            args.feature_counts, sep="\t", comment="#", dtype=str, keep_default_na=False
        )
        required_count_columns = {"Geneid", args.expression_sample_column}
        if not required_count_columns.issubset(counts.columns):
            raise ValueError(
                "FeatureCounts table lacks required columns: "
                f"{sorted(required_count_columns - set(counts.columns))}"
            )
        if counts["Geneid"].duplicated().any():
            raise ValueError("FeatureCounts table contains duplicate Geneid values")
        counts = counts.set_index("Geneid")
        missing_expression = reference - set(counts.index)
        if missing_expression:
            raise ValueError(
                f"FeatureCounts table lacks {len(missing_expression)} reference genes"
            )
        values = pd.to_numeric(
            counts.loc[analysis_order, args.expression_sample_column], errors="raise"
        )
        if (values < 0).any() or (values % 1 != 0).any():
            raise ValueError("Selected expression values must be non-negative integer raw counts")
        expression = pd.DataFrame(
            {
                "reference_gene_id": analysis_order,
                "leaf_raw_count": values.astype("int64").to_numpy(),
            }
        )

        cluster_members = read_cluster_members(args.clusters)
        unknown_members = cluster_members - reference
        if unknown_members:
            raise ValueError(
                f"CD-HIT membership contains {len(unknown_members)} non-reference IDs"
            )
        missing_cluster = sorted(analysis_reference - cluster_members)

        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(
                prefix=f".{args.output_dir.name}.tmp.", dir=args.output_dir.parent
            )
        )
        try:
            outputs = {
                "resolved_nonshared_unit_loss_table.tsv": resolved,
                "nonshared_resolved_denominator_audit.tsv": audit,
                "nonshared_reference_leaf_raw_counts.tsv": expression,
            }
            for name, frame in outputs.items():
                write_tsv(staging / name, frame)
            (staging / "nonshared_cdhit_missing_reference_ids.txt").write_text(
                "".join(f"{gene}\n" for gene in missing_cluster), encoding="utf-8"
            )
            output_names = [*outputs, "nonshared_cdhit_missing_reference_ids.txt"]
            checksum_rows = []
            for name in output_names:
                path = staging / name
                checksum_rows.append(
                    {"file": name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                )
            with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["file", "bytes", "sha256"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(checksum_rows)
            manifest = {
                "schema_version": 1,
                "status": "PASS",
                "workflow": "primary_resolved_downstream_inputs",
                "denominator_policy": (
                    "shared positive-complete reference genes excluded; "
                    "deleted plus strictly supported pseudogenized plus retained; uncertain excluded"
                ),
                "expression_measurement": "raw_count",
                "expression_tissue": "leaf",
                "expression_sample_column": args.expression_sample_column,
                "reference_gene_count": len(reference),
                "shared_positive_gene_count": len(shared_positive),
                "nonshared_analysis_gene_count": len(analysis_reference),
                "assembly_unit_count": len(unit_ids),
                "complete_loss_rows": len(losses),
                "resolved_loss_rows": len(resolved),
                "uncertain_loss_rows": int(
                    (nonshared_losses["classification"] == "uncertain").sum()
                ),
                "expression_gene_count": len(expression),
                "cdhit_clustered_reference_genes": len(cluster_members),
                "cdhit_clustered_nonshared_analysis_genes": len(
                    analysis_reference & cluster_members
                ),
                "cdhit_missing_reference_genes": len(missing_cluster),
                "inputs": [
                    {
                        "label": label,
                        "basename": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                    for label, path in zip(
                        [
                            "loss_matrix", "unit_ledger", "reference_cds",
                            "shared_positive_genes",
                            "feature_counts", "cdhit_clusters",
                        ],
                        inputs,
                    )
                ],
                "outputs": checksum_rows,
            }
            (staging / "run_manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            os.replace(staging, args.output_dir)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        f"PASS\t{args.output_dir}\tunits={args.expected_units}\t"
        f"reference_genes={args.expected_reference_genes}\t"
        f"cdhit_missing={len(missing_cluster)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
