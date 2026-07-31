#!/usr/bin/env python3
"""Build one old/new complete loss matrix from the same local evidence rules.

SynOrths anchors are retained.  Every non-anchor candidate is taken from the
same bilateral-anchor candidate builder and the same Miniprot search.  A
callable locus without a qualifying local alignment is deleted.  A local
alignment is pseudogenized only when it contains an explicit frameshift or
in-frame stop *and* passes the stricter disruption-quality gate.  All other
cases are uncertain.  The resulting table retains exact position semantics so
that positive loss versus retained can be modelled on a resolved denominator.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


class MatrixError(RuntimeError):
    pass


CONFIG_COLUMNS = (
    "unit",
    "source_group",
    "synorth_pairs",
    "synorth_reference_column_1based",
    "candidate_dir",
    "uniform_output_dir",
)
OUTPUT_COLUMNS = (
    "reference_gene_id",
    "assembly_unit_id",
    "source_group",
    "classification",
    "callable",
    "positive_loss",
    "evidence_source",
    "evidence_reason",
    "target_chromosome",
    "position_start_1based",
    "position_end_1based",
    "position_midpoint_1based",
    "coordinate_semantics",
    "resolved_for_spatial_model",
    "callable_interval_chromosome",
    "callable_interval_start_1based",
    "callable_interval_end_1based",
    "callable_interval_midpoint_1based",
    "query_coverage",
    "exact_alignment_identity",
    "alignment_score",
    "frameshift_events",
    "inframe_stop_codons",
    "disruption_supported",
)
DISRUPTION_COLUMNS = (
    "reference_gene_id",
    "assembly_unit_id",
    "source_group",
    "target_chromosome",
    "position_start_1based",
    "position_end_1based",
    "query_coverage",
    "exact_alignment_identity",
    "alignment_score",
    "frameshift_events",
    "inframe_stop_codons",
    "evidence_reason",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.resolve().open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (source.stat().st_size == 0 and not allow_empty):
        raise MatrixError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MatrixError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise MatrixError(f"{path}: JSON root is not an object")
    return value


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise MatrixError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise MatrixError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size == 0):
        raise MatrixError(f"missing input file: {value!r}")
    return path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not reader.fieldnames:
        raise MatrixError(f"{path.name}: missing header")
    return list(reader.fieldnames), rows


def read_config(path: Path) -> list[dict[str, str]]:
    columns, rows = read_tsv(path)
    if tuple(columns) != CONFIG_COLUMNS or not rows:
        raise MatrixError("config columns/rows differ from exact schema")
    units = [row["unit"] for row in rows]
    if len(units) != len(set(units)) or any(not row[key] for row in rows for key in CONFIG_COLUMNS):
        raise MatrixError("config has duplicate units or empty values")
    if any(row["source_group"] not in {"legacy", "new"} for row in rows):
        raise MatrixError("source_group must be legacy or new")
    return rows


def read_fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            identifier = raw[1:].strip().split()[0]
            if not identifier:
                raise MatrixError(f"{path.name}:{line_number}: empty FASTA ID")
            identifiers.append(identifier)
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise MatrixError(f"{path.name}: no or duplicate FASTA IDs")
    return identifiers


def read_synorths(path: Path, reference_column: int, reference_ids: set[str]) -> dict[str, list[tuple[str, int, int]]]:
    if reference_column not in {1, 5}:
        raise MatrixError("SynOrths reference column must be 1 or 5")
    target_offset = 4 if reference_column == 1 else 0
    anchors: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise MatrixError(f"{path.name}:{line_number}: fewer than 8 columns")
            gene = fields[reference_column - 1]
            if gene not in reference_ids:
                raise MatrixError(f"{path.name}:{line_number}: unknown reference gene {gene!r}")
            chromosome = fields[target_offset + 1]
            try:
                start = int(fields[target_offset + 2])
                end = int(fields[target_offset + 3])
            except ValueError as error:
                raise MatrixError(f"{path.name}:{line_number}: invalid target coordinates") from error
            start, end = min(start, end), max(start, end)
            if not chromosome or start < 1:
                raise MatrixError(f"{path.name}:{line_number}: invalid target locus")
            anchors[gene].append((chromosome, start, end))
    if not anchors:
        raise MatrixError(f"{path.name}: no SynOrths anchors")
    return dict(anchors)


def validate_bundle(directory: Path, unit: str, table_name: str) -> tuple[Path, dict[str, object]]:
    manifest_path = directory / "run_manifest.json"
    table_path = directory / table_name
    checksum_path = directory / "checksums.tsv"
    for path in (manifest_path, table_path, checksum_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise MatrixError(f"{unit}: incomplete bundle {directory}")
    manifest = strict_json(manifest_path)
    if manifest.get("status") != "PASS" or manifest.get("unit") != unit:
        raise MatrixError(f"{unit}: bundle is not exact PASS")
    with checksum_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    match = [row for row in rows if row.get("file") == table_name]
    if len(match) != 1:
        raise MatrixError(f"{unit}: checksum row missing for {table_name}")
    observed = binding(table_path, allow_empty=True)
    if str(observed["bytes"]) != match[0].get("bytes") or observed["sha256"] != match[0].get("sha256"):
        raise MatrixError(f"{unit}: checksum mismatch for {table_name}")
    return table_path, manifest


def read_candidates(path: Path, unit: str) -> dict[str, dict[str, str]]:
    columns, rows = read_tsv(path)
    required = {
        "unit", "reference_gene", "callable", "callability_reason", "target_chromosome",
        "target_interval_start_1based", "target_interval_end_1based",
    }
    if not rows or not required.issubset(columns):
        raise MatrixError(f"{unit}: candidate schema/rows invalid")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        gene = row["reference_gene"]
        if row["unit"] != unit or not gene or gene in result or row["callable"] not in {"true", "false"}:
            raise MatrixError(f"{unit}: candidate unit/gene/boolean closure failed")
        result[gene] = row
    return result


def read_states(path: Path, unit: str) -> dict[str, dict[str, str]]:
    columns, rows = read_tsv(path)
    required = {
        "unit", "reference_gene", "callable", "classification", "qualifying_local_alignment",
        "alignment_target_start_1based", "alignment_target_end_1based", "query_coverage",
        "exact_alignment_identity", "alignment_score", "frameshift_events", "inframe_stop_codons",
    }
    if not rows or not required.issubset(columns):
        raise MatrixError(f"{unit}: uniform state schema/rows invalid")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        gene = row["reference_gene"]
        if row["unit"] != unit or not gene or gene in result:
            raise MatrixError(f"{unit}: uniform state unit/gene closure failed")
        result[gene] = row
    return result


def numeric(value: str, context: str, kind=float):
    try:
        return kind(value)
    except ValueError as error:
        raise MatrixError(f"{context}: invalid number {value!r}") from error


def retained_row(unit: str, source_group: str, gene: str, loci: list[tuple[str, int, int]]) -> dict[str, str]:
    unique = sorted(set(loci))
    if len(unique) == 1:
        chromosome, start, end = unique[0]
        position = (start + end) / 2
        resolved, semantics = "true", "exact_synorth_target_locus_midpoint"
    else:
        chromosome, start, end, position = "", "", "", ""
        resolved, semantics = "false", "multiple_synorth_target_loci_no_unique_position"
    return {
        "reference_gene_id": gene,
        "assembly_unit_id": unit,
        "source_group": source_group,
        "classification": "retained",
        "callable": "true",
        "positive_loss": "false",
        "evidence_source": "exact_synorth_anchor",
        "evidence_reason": "exact_synorth_anchor_present",
        "target_chromosome": chromosome,
        "position_start_1based": str(start),
        "position_end_1based": str(end),
        "position_midpoint_1based": f"{position:.12g}" if position != "" else "",
        "coordinate_semantics": semantics,
        "resolved_for_spatial_model": resolved,
        "callable_interval_chromosome": "",
        "callable_interval_start_1based": "",
        "callable_interval_end_1based": "",
        "callable_interval_midpoint_1based": "",
        "query_coverage": "",
        "exact_alignment_identity": "",
        "alignment_score": "",
        "frameshift_events": "",
        "inframe_stop_codons": "",
        "disruption_supported": "false",
    }


def candidate_row(
    unit: str,
    source_group: str,
    gene: str,
    candidate: dict[str, str],
    state: dict[str, str],
    *,
    disruption_coverage: float,
    disruption_identity: float,
    disruption_score: int,
) -> dict[str, str]:
    if candidate["callable"] != state["callable"]:
        raise MatrixError(f"{unit}/{gene}: candidate/state callability differs")
    callable_value = candidate["callable"] == "true"
    raw_class = state["classification"]
    disruption = False
    if raw_class == "deleted":
        classification, reason = "deleted", "no_qualifying_local_protein_genome_alignment"
    elif raw_class == "pseudogenized":
        values = (
            numeric(state["query_coverage"], f"{unit}/{gene}:coverage"),
            numeric(state["exact_alignment_identity"], f"{unit}/{gene}:identity"),
            numeric(state["alignment_score"], f"{unit}/{gene}:score", int),
            numeric(state["frameshift_events"], f"{unit}/{gene}:frameshifts", int),
            numeric(state["inframe_stop_codons"], f"{unit}/{gene}:stops", int),
        )
        disruption = (
            values[0] >= disruption_coverage
            and values[1] >= disruption_identity
            and values[2] >= disruption_score
            and (values[3] > 0 or values[4] > 0)
        )
        if disruption:
            classification, reason = "pseudogenized", "high_confidence_frameshift_or_inframe_stop"
        else:
            classification, reason = "uncertain", "disruptive_tag_below_strict_quality_gate"
    elif raw_class == "uncertain":
        classification = "uncertain"
        reason = state.get("evidence_reason", "local_sequence_without_disruptive_event")
    else:
        raise MatrixError(f"{unit}/{gene}: unsupported raw class {raw_class!r}")
    if classification in {"deleted", "pseudogenized"} and not callable_value:
        raise MatrixError(f"{unit}/{gene}: positive call is non-callable")

    if classification == "deleted":
        chromosome = candidate["target_chromosome"]
        start = candidate["target_interval_start_1based"]
        end = candidate["target_interval_end_1based"]
        semantics = "midpoint_of_bilateral_synorth_bounded_interval"
    elif classification == "pseudogenized":
        chromosome = candidate["target_chromosome"]
        start = state["alignment_target_start_1based"]
        end = state["alignment_target_end_1based"]
        semantics = "midpoint_of_disrupted_local_miniprot_alignment"
    else:
        chromosome = start = end = ""
        semantics = "excluded_uncertain"
    if start and end:
        start_number = numeric(start, f"{unit}/{gene}:position start", int)
        end_number = numeric(end, f"{unit}/{gene}:position end", int)
        midpoint = f"{(start_number + end_number) / 2:.12g}"
    else:
        midpoint = ""
    interval_chromosome = candidate["target_chromosome"] if callable_value else ""
    interval_start = candidate["target_interval_start_1based"] if callable_value else ""
    interval_end = candidate["target_interval_end_1based"] if callable_value else ""
    if interval_start and interval_end:
        interval_midpoint = f"{(int(interval_start) + int(interval_end)) / 2:.12g}"
    else:
        interval_midpoint = ""
    return {
        "reference_gene_id": gene,
        "assembly_unit_id": unit,
        "source_group": source_group,
        "classification": classification,
        "callable": str(callable_value).lower(),
        "positive_loss": str(classification in {"deleted", "pseudogenized"}).lower(),
        "evidence_source": "uniform_bilateral_anchor_miniprot",
        "evidence_reason": reason,
        "target_chromosome": chromosome,
        "position_start_1based": start,
        "position_end_1based": end,
        "position_midpoint_1based": midpoint,
        "coordinate_semantics": semantics,
        "resolved_for_spatial_model": str(classification in {"retained", "deleted", "pseudogenized"}).lower(),
        "callable_interval_chromosome": interval_chromosome,
        "callable_interval_start_1based": interval_start,
        "callable_interval_end_1based": interval_end,
        "callable_interval_midpoint_1based": interval_midpoint,
        "query_coverage": state["query_coverage"],
        "exact_alignment_identity": state["exact_alignment_identity"],
        "alignment_score": state["alignment_score"],
        "frameshift_events": state["frameshift_events"],
        "inframe_stop_codons": state["inframe_stop_codons"],
        "disruption_supported": str(disruption).lower(),
    }


def unscoped_row(unit: str, source_group: str, gene: str) -> dict[str, str]:
    row = {column: "" for column in OUTPUT_COLUMNS}
    row.update(
        {
            "reference_gene_id": gene,
            "assembly_unit_id": unit,
            "source_group": source_group,
            "classification": "uncertain",
            "callable": "false",
            "positive_loss": "false",
            "evidence_source": "outside_bilateral_candidate_scope",
            "evidence_reason": "no_exact_anchor_or_candidate_interval",
            "coordinate_semantics": "excluded_uncertain",
            "resolved_for_spatial_model": "false",
            "disruption_supported": "false",
        }
    )
    return row


def write_tsv(path: Path, rows, columns) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-disruption-query-coverage", type=float, default=0.80)
    parser.add_argument("--minimum-disruption-identity", type=float, default=0.70)
    parser.add_argument("--minimum-disruption-alignment-score", type=int, default=100)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        if not 0 < args.minimum_disruption_query_coverage <= 1 or not 0 < args.minimum_disruption_identity <= 1:
            raise MatrixError("invalid disruption fraction threshold")
        if args.minimum_disruption_alignment_score < 0:
            raise MatrixError("invalid disruption score threshold")
        root = args.data_root.resolve()
        config_rows = read_config(args.config)
        reference_protein = args.reference_protein.resolve()
        reference_genes = read_fasta_ids(reference_protein)
        reference_set = set(reference_genes)
        output = args.output_dir.resolve()
        if output.exists():
            raise MatrixError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))

        matrix_path = staging / "complete_unit_loss_matrix.tsv"
        summaries: list[dict[str, object]] = []
        input_audit: list[dict[str, object]] = []
        disruption_path = staging / "strict_disruption_calls.tsv"
        with (
            matrix_path.open("w", encoding="utf-8", newline="") as handle,
            disruption_path.open("w", encoding="utf-8", newline="") as disruption_handle,
        ):
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
            disruption_writer = csv.DictWriter(
                disruption_handle, fieldnames=DISRUPTION_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            disruption_writer.writeheader()
            for config in config_rows:
                unit = config["unit"]
                source_group = config["source_group"]
                synorth_path = resolve(root, config["synorth_pairs"], file=True)
                candidate_dir = resolve(root, config["candidate_dir"])
                uniform_dir = resolve(root, config["uniform_output_dir"])
                candidate_path, candidate_manifest = validate_bundle(candidate_dir, unit, "candidates.tsv")
                state_path, state_manifest = validate_bundle(uniform_dir, unit, "uniform_candidate_loss_states.tsv")
                anchors = read_synorths(
                    synorth_path, int(config["synorth_reference_column_1based"]), reference_set
                )
                candidates = read_candidates(candidate_path, unit)
                states = read_states(state_path, unit)
                if set(candidates) != set(states) or set(anchors).intersection(candidates):
                    raise MatrixError(f"{unit}: anchor/candidate/state disjoint closure failed")
                outside = reference_set.difference(anchors).difference(candidates)
                counts: Counter[str] = Counter()
                resolved = 0
                disruption_count = 0
                for gene in reference_genes:
                    if gene in anchors:
                        result = retained_row(unit, source_group, gene, anchors[gene])
                    elif gene in candidates:
                        result = candidate_row(
                            unit,
                            source_group,
                            gene,
                            candidates[gene],
                            states[gene],
                            disruption_coverage=args.minimum_disruption_query_coverage,
                            disruption_identity=args.minimum_disruption_identity,
                            disruption_score=args.minimum_disruption_alignment_score,
                        )
                    else:
                        result = unscoped_row(unit, source_group, gene)
                    writer.writerow(result)
                    if result["classification"] == "pseudogenized":
                        disruption_writer.writerow({column: result[column] for column in DISRUPTION_COLUMNS})
                    counts[result["classification"]] += 1
                    resolved += result["resolved_for_spatial_model"] == "true"
                    disruption_count += result["disruption_supported"] == "true"
                summaries.append(
                    {
                        "assembly_unit_id": unit,
                        "source_group": source_group,
                        "reference_gene_count": len(reference_genes),
                        "retained": counts["retained"],
                        "deleted": counts["deleted"],
                        "pseudogenized": counts["pseudogenized"],
                        "uncertain": counts["uncertain"],
                        "positive_loss": counts["deleted"] + counts["pseudogenized"],
                        "strict_disruption_supported": disruption_count,
                        "resolved_for_spatial_model": resolved,
                        "outside_candidate_scope": len(outside),
                    }
                )
                input_audit.append(
                    {
                        "unit": unit,
                        "synorth_pairs": binding(synorth_path),
                        "candidate_manifest": binding(candidate_dir / "run_manifest.json"),
                        "uniform_manifest": binding(uniform_dir / "run_manifest.json"),
                        "candidate_metrics": candidate_manifest.get("metrics"),
                        "uniform_metrics": state_manifest.get("metrics"),
                    }
                )
        summary_columns = tuple(summaries[0])
        write_tsv(staging / "unit_summary.tsv", summaries, summary_columns)
        definitions = [
            {"classification": "retained", "positive_loss": "false", "definition": "exact SynOrths anchor present"},
            {"classification": "deleted", "positive_loss": "true", "definition": "callable bilateral interval and no qualifying local Miniprot alignment"},
            {"classification": "pseudogenized", "positive_loss": "true", "definition": "callable local alignment with frameshift or in-frame stop and strict alignment-quality support"},
            {"classification": "uncertain", "positive_loss": "false", "definition": "non-callable, local sequence without a supported disruptive event, or outside candidate scope"},
        ]
        write_tsv(staging / "classification_definitions.tsv", definitions, ("classification", "positive_loss", "definition"))
        manifest = {
            "schema_version": 1,
            "workflow": "uniform_old_new_complete_loss_matrix",
            "status": "PASS",
            "finished_at_utc": utc_now(),
            "parameters": {
                "reference_gene_count": len(reference_genes),
                "assembly_unit_count": len(config_rows),
                "matrix_row_count": len(reference_genes) * len(config_rows),
                "local_sequence_gate": {"query_coverage": 0.50, "exact_identity": 0.50, "alignment_score": 50},
                "strict_disruption_gate": {
                    "query_coverage": args.minimum_disruption_query_coverage,
                    "exact_identity": args.minimum_disruption_identity,
                    "alignment_score": args.minimum_disruption_alignment_score,
                    "required_event": "Miniprot fs>0 or st>0",
                },
                "positive_classes": ["deleted", "pseudogenized"],
                "uncertain_is_never_positive": True,
                "raw_read_validation": "not_available; disruptions are assembly-sequence-supported",
            },
            "inputs": {
                "config": binding(args.config),
                "reference_protein": binding(reference_protein),
                "unit_bundles": input_audit,
            },
            "metrics": {"unit_summaries": summaries},
        }
        (staging / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(json.dumps({"status": "PASS", "rows": len(reference_genes) * len(config_rows), "output": str(output)}, sort_keys=True))
        return 0
    except (OSError, csv.Error, ValueError, MatrixError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    finally:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
