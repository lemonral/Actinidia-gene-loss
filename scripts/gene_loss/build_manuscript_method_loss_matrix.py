#!/usr/bin/env python3
"""Build a 23-unit manuscript-method loss matrix plus conservative refinements.

The first classification layer reproduces the manuscript rule without changing
its interpretation: an exact SynOrths anchor is retained; a missing-gene
candidate with a qualifying genome-wide tBLASTX hit is decayed; and a candidate
without such a hit is deleted.  Rows outside the historical candidate scope are
``not_called_loss`` rather than being forced into a loss class.

The second layer never changes the manuscript class.  It partitions ``decayed``
rows using already-generated Miniprot evidence.  Explicit frameshift/stop calls
and alignment-end truncation candidates are reported separately; all remaining
rows retain an explicit unresolved label.  No new sequence search is run.
"""

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


class BuildError(RuntimeError):
    """Raised when an input or evidence closure check fails."""


UNIFORM_REQUIRED = {
    "reference_gene_id",
    "assembly_unit_id",
    "source_group",
    "classification",
    "callable",
    "evidence_reason",
    "query_coverage",
    "exact_alignment_identity",
    "alignment_score",
    "frameshift_events",
    "inframe_stop_codons",
    "disruption_supported",
}

MINIPROT_REQUIRED = {
    "unit",
    "reference_gene",
    "qualifying_local_alignment",
    "query_length_aa",
    "query_aligned_start_0based",
    "query_aligned_end_0based_exclusive",
    "query_coverage",
    "frameshift_events",
    "inframe_stop_codons",
}

OUTPUT_FIELDS = (
    "reference_gene_id",
    "assembly_unit_id",
    "source_group",
    "manuscript_classification",
    "manuscript_positive_loss",
    "manuscript_rule",
    "callable",
    "uniform_classification",
    "uniform_evidence_reason",
    "refined_decayed_cause",
    "refined_cause_evidence_level",
    "query_coverage",
    "exact_alignment_identity",
    "alignment_score",
    "frameshift_events",
    "inframe_stop_codons",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise BuildError(f"missing or empty {label}: {path}")
    return path


def rows(path: Path, required: set[str] | None = None):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BuildError(f"missing TSV header: {path}")
        if required:
            missing = required.difference(reader.fieldnames)
            if missing:
                raise BuildError(f"{path} missing columns: {sorted(missing)}")
        yield from reader


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def load_legacy(path: Path) -> dict[tuple[str, str], str]:
    required = {"reference_gene_id", "assembly_unit_id", "classification"}
    result: dict[tuple[str, str], str] = {}
    translation = {
        "retained": "retained",
        "pseudogenized": "decayed",
        "deleted": "deleted",
        "not_called_loss": "not_called_loss",
    }
    for row in rows(path, required):
        key = (row["assembly_unit_id"], row["reference_gene_id"])
        value = translation.get(row["classification"])
        if value is None or key in result:
            raise BuildError(f"invalid or duplicate legacy row: {key}")
        result[key] = value
    return result


def load_new_states(config: Path, root: Path) -> tuple[dict[tuple[str, str], str], list[Path]]:
    result: dict[tuple[str, str], str] = {}
    inputs: list[Path] = []
    for config_row in rows(config, {"unit", "output_dir"}):
        unit = config_row["unit"]
        path = require_file(resolve(root, config_row["output_dir"]) / "loss_states.tsv", "new loss states")
        inputs.append(path)
        for row in rows(path, {"unit", "reference_gene", "historical_reproduction_state"}):
            if row["unit"] != unit:
                raise BuildError(f"unit mismatch in {path}")
            key = (unit, row["reference_gene"])
            value = row["historical_reproduction_state"]
            if value not in {"decayed", "deleted"} or key in result:
                raise BuildError(f"invalid or duplicate new historical row: {key}")
            result[key] = value
    return result, inputs


def load_miniprot_states(config: Path, root: Path) -> tuple[dict[tuple[str, str], dict[str, str]], list[Path]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    inputs: list[Path] = []
    for config_row in rows(config, {"unit", "uniform_output_dir"}):
        unit = config_row["unit"]
        path = require_file(
            resolve(root, config_row["uniform_output_dir"]) / "uniform_candidate_loss_states.tsv",
            "Miniprot candidate states",
        )
        inputs.append(path)
        for row in rows(path, MINIPROT_REQUIRED):
            if row["unit"] != unit:
                raise BuildError(f"unit mismatch in {path}")
            key = (unit, row["reference_gene"])
            if key in result:
                raise BuildError(f"duplicate Miniprot state: {key}")
            result[key] = row
    return result, inputs


def integer(value: str) -> int:
    return 0 if value == "" else int(value)


def refined_decayed_cause(
    manuscript_classification: str,
    uniform: dict[str, str],
    miniprot: dict[str, str] | None,
    *,
    terminal_missing_fraction: float,
) -> tuple[str, str]:
    """Return a cause label and its evidence level without changing old class."""
    if manuscript_classification == "retained":
        return "not_applicable_retained", "exact_synorth"
    if manuscript_classification == "not_called_loss":
        return "not_called_outside_historical_scope", "not_callable_as_loss"

    # Keep this orthogonal to the old decayed/deleted split.  Miniprot can find
    # a protein-to-genome alignment missed by the old nucleotide tBLASTX gate;
    # that extra evidence is worth retaining without rewriting the old class.
    frameshifts = integer(uniform["frameshift_events"])
    stops = integer(uniform["inframe_stop_codons"])
    strict = uniform["disruption_supported"] == "true"
    if strict and frameshifts and stops:
        return "frameshift_and_stop_supported", "explicit_coding_disruption"
    if strict and frameshifts:
        return "frameshift_supported", "explicit_coding_disruption"
    if strict and stops:
        return "stop_supported", "explicit_coding_disruption"
    if frameshifts or stops:
        return "frameshift_or_stop_below_strict_quality_gate", "candidate_only"

    if manuscript_classification == "deleted":
        return "no_qualifying_genomewide_tblastx_hit", "manuscript_threshold"

    if miniprot is not None and miniprot["qualifying_local_alignment"] == "true":
        qlen = integer(miniprot["query_length_aa"])
        qstart = integer(miniprot["query_aligned_start_0based"])
        qend = integer(miniprot["query_aligned_end_0based_exclusive"])
        if qlen <= 0 or not 0 <= qstart <= qend <= qlen:
            raise BuildError("invalid Miniprot query interval")
        n_missing = qstart / qlen
        c_missing = (qlen - qend) / qlen
        n_flag = n_missing >= terminal_missing_fraction
        c_flag = c_missing >= terminal_missing_fraction
        if n_flag and c_flag:
            return "both_terminal_alignment_truncation_candidate", "alignment_candidate_only"
        if n_flag:
            return "n_terminal_alignment_truncation_candidate", "alignment_candidate_only"
        if c_flag:
            return "c_terminal_alignment_truncation_candidate", "alignment_candidate_only"
        if float(miniprot["query_coverage"]) < 0.80:
            return "partial_local_alignment_other_candidate", "alignment_candidate_only"
        return "local_sequence_no_explicit_coding_disruption", "sequence_detected_only"

    if uniform["callable"] != "true":
        return "genomewide_tblastx_hit_noncallable_local_locus", "sequence_detected_only"
    return "genomewide_tblastx_hit_without_local_miniprot_support", "sequence_detected_only"


def manuscript_classification(
    row: dict[str, str],
    legacy: dict[tuple[str, str], str],
    new: dict[tuple[str, str], str],
) -> str:
    key = (row["assembly_unit_id"], row["reference_gene_id"])
    if row["source_group"] == "legacy":
        if key not in legacy:
            raise BuildError(f"legacy historical matrix is incomplete: {key}")
        return legacy[key]
    if row["source_group"] != "new":
        raise BuildError(f"unknown source_group: {row['source_group']}")
    if row["classification"] == "retained":
        if key in new:
            raise BuildError(f"new retained gene also occurs in candidate states: {key}")
        return "retained"
    return new.get(key, "not_called_loss")


def write_tsv(path: Path, output_rows: list[dict[str, object]], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)


def build(args: argparse.Namespace) -> None:
    root = args.data_root.resolve()
    uniform_path = require_file(args.uniform_matrix.resolve(), "uniform matrix")
    legacy_path = require_file(args.legacy_historical_matrix.resolve(), "legacy historical matrix")
    legacy = load_legacy(legacy_path)
    new, new_inputs = load_new_states(args.new_search_config.resolve(), root)
    miniprot, miniprot_inputs = load_miniprot_states(args.miniprot_config.resolve(), root)

    output = args.output_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    manuscript_counts: Counter[tuple[str, str]] = Counter()
    cause_counts: Counter[tuple[str, str]] = Counter()
    crosswalk_counts: Counter[tuple[str, str, str]] = Counter()
    unit_genes: Counter[str] = Counter()
    seen: set[tuple[str, str]] = set()
    matrix_path = staging / "manuscript_method_unit_gene_matrix.tsv.gz"
    try:
        with gzip.open(matrix_path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=OUTPUT_FIELDS, lineterminator="\n")
            writer.writeheader()
            for row in rows(uniform_path, UNIFORM_REQUIRED):
                key = (row["assembly_unit_id"], row["reference_gene_id"])
                if key in seen:
                    raise BuildError(f"duplicate uniform matrix row: {key}")
                seen.add(key)
                old = manuscript_classification(row, legacy, new)
                cause, level = refined_decayed_cause(
                    old, row, miniprot.get(key), terminal_missing_fraction=args.terminal_missing_fraction
                )
                manuscript_counts[(row["assembly_unit_id"], old)] += 1
                cause_counts[(row["assembly_unit_id"], cause)] += 1
                crosswalk_counts[(old, row["classification"], cause)] += 1
                unit_genes[row["assembly_unit_id"]] += 1
                writer.writerow(
                    {
                        "reference_gene_id": row["reference_gene_id"],
                        "assembly_unit_id": row["assembly_unit_id"],
                        "source_group": row["source_group"],
                        "manuscript_classification": old,
                        "manuscript_positive_loss": str(old in {"decayed", "deleted"}).lower(),
                        "manuscript_rule": "synorth_anchor_else_genomewide_tblastx_50pct_bitscore50_evalue_lt_1e-5",
                        "callable": row["callable"],
                        "uniform_classification": row["classification"],
                        "uniform_evidence_reason": row["evidence_reason"],
                        "refined_decayed_cause": cause,
                        "refined_cause_evidence_level": level,
                        "query_coverage": row["query_coverage"],
                        "exact_alignment_identity": row["exact_alignment_identity"],
                        "alignment_score": row["alignment_score"],
                        "frameshift_events": row["frameshift_events"],
                        "inframe_stop_codons": row["inframe_stop_codons"],
                    }
                )

        if len(unit_genes) != args.expected_units or any(value != args.expected_genes for value in unit_genes.values()):
            raise BuildError(f"unexpected unit/gene grid: {dict(unit_genes)}")

        summary_rows: list[dict[str, object]] = []
        for unit in sorted(unit_genes):
            counts = {name: manuscript_counts[(unit, name)] for name in ("retained", "decayed", "deleted", "not_called_loss")}
            denominator = counts["retained"] + counts["decayed"] + counts["deleted"]
            summary_rows.append(
                {
                    "assembly_unit_id": unit,
                    **counts,
                    "manuscript_positive_loss": counts["decayed"] + counts["deleted"],
                    "manuscript_resolved_denominator": denominator,
                    "manuscript_loss_rate": f"{(counts['decayed'] + counts['deleted']) / denominator:.12g}",
                }
            )
        write_tsv(
            staging / "unit_summary.tsv",
            summary_rows,
            (
                "assembly_unit_id", "retained", "decayed", "deleted", "not_called_loss",
                "manuscript_positive_loss", "manuscript_resolved_denominator", "manuscript_loss_rate",
            ),
        )

        cause_rows = [
            {"assembly_unit_id": unit, "refined_decayed_cause": cause, "row_count": count}
            for (unit, cause), count in sorted(cause_counts.items())
        ]
        write_tsv(staging / "refined_cause_summary.tsv", cause_rows, ("assembly_unit_id", "refined_decayed_cause", "row_count"))
        crosswalk_rows = [
            {
                "manuscript_classification": old,
                "uniform_classification": uniform,
                "refined_decayed_cause": cause,
                "row_count": count,
            }
            for (old, uniform, cause), count in sorted(crosswalk_counts.items())
        ]
        write_tsv(
            staging / "classification_crosswalk.tsv",
            crosswalk_rows,
            ("manuscript_classification", "uniform_classification", "refined_decayed_cause", "row_count"),
        )
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_MANUSCRIPT_METHOD_WITH_CONSERVATIVE_REFINEMENT",
            "unit_count": len(unit_genes),
            "reference_gene_count": args.expected_genes,
            "matrix_rows": sum(unit_genes.values()),
            "manuscript_rule": {
                "retained": "exact SynOrths anchor",
                "decayed": "candidate has genome-wide tBLASTX hit with identity>=50%, bitscore>=50, evalue<1e-5, no length minimum",
                "deleted": "candidate lacks such a qualifying genome-wide tBLASTX hit",
                "not_called_loss": "outside the historical missing-gene candidate scope",
            },
            "refinement_policy": {
                "does_not_change_manuscript_classification": True,
                "terminal_alignment_truncation_candidate_threshold": args.terminal_missing_fraction,
                "explicit_disruption_source": "existing Miniprot fs/st tags",
                "unsupported_mechanisms_not_inferred": [
                    "start_codon_loss", "terminal_stop_loss", "splice_site_disruption",
                    "exon_loss", "gene_fission", "gene_fusion", "TE_insertion", "epigenetic_silencing",
                ],
            },
            "inputs": {
                "uniform_matrix": binding(uniform_path),
                "legacy_historical_matrix": binding(legacy_path),
                "new_loss_state_files": [binding(path) for path in new_inputs],
                "miniprot_state_files": [binding(path) for path in miniprot_inputs],
            },
        }
        outputs = [matrix_path, staging / "unit_summary.tsv", staging / "refined_cause_summary.tsv", staging / "classification_crosswalk.tsv"]
        manifest["outputs"] = [binding(path) for path in outputs]
        (staging / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        if output.exists():
            raise BuildError(f"refusing to overwrite existing output: {output}")
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--uniform-matrix", type=Path, required=True)
    parser.add_argument("--legacy-historical-matrix", type=Path, required=True)
    parser.add_argument("--new-search-config", type=Path, required=True)
    parser.add_argument("--miniprot-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--terminal-missing-fraction", type=float, default=0.20)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-genes", type=int, default=35547)
    args = parser.parse_args()
    if not 0 < args.terminal_missing_fraction < 0.5:
        parser.error("--terminal-missing-fraction must be between 0 and 0.5")
    return args


if __name__ == "__main__":
    build(parse_args())
