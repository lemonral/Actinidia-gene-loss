"""tBLASTX schema validation and pseudogene/deletion classification.

The archived repository contains both a six-column parser and existing
12-column BLAST outfmt-6 files.  This module never silently applies one parser
to the other.  It records the detected schema and writes the best qualifying
hit used for every classification decision.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .gff import iter_fasta
from .io_utils import SchemaError, format_number, parse_float, parse_int, read_id_file, read_tsv, write_tsv


STANDARD_12 = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
HEADER_ALIASES = {
    "query_id": "qseqid", "subject_id": "sseqid", "identity": "pident",
    "alignment_length": "length", "bit_score": "bitscore", "score": "bitscore",
    "e_value": "evalue", "query_start": "qstart", "query_end": "qend",
    "subject_start": "sstart", "subject_end": "send",
}


@dataclass(frozen=True)
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    evalue: float
    bitscore: float
    qstart: int | None = None
    qend: int | None = None
    sstart: int | None = None
    send: int | None = None
    raw_line_number: int = 0

    @property
    def subject_start(self) -> int | None:
        if self.sstart is None or self.send is None:
            return None
        return min(self.sstart, self.send)

    @property
    def subject_end(self) -> int | None:
        if self.sstart is None or self.send is None:
            return None
        return max(self.sstart, self.send)


def _split(raw: str) -> list[str]:
    return raw.rstrip("\n").split("\t") if "\t" in raw else raw.split()


def _looks_like_header(values: list[str]) -> bool:
    lowered = {value.strip().lower() for value in values}
    return bool(lowered & {"qseqid", "query_id", "sseqid", "subject_id", "pident", "identity", "evalue", "bitscore"})


def _canonical_header(values: list[str], source: str) -> dict[str, int]:
    canonical: dict[str, int] = {}
    for index, value in enumerate(values):
        name = HEADER_ALIASES.get(value.strip().lower(), value.strip().lower())
        canonical[name] = index
    required = {"qseqid", "sseqid", "pident", "length", "evalue", "bitscore"}
    missing = sorted(required - set(canonical))
    if missing:
        raise SchemaError(f"{source}: headered BLAST table is missing required columns: {', '.join(missing)}")
    return canonical


def _six_column_mapping(values: list[str], source: str, requested: str) -> tuple[int, int]:
    """Return (bitscore_index, evalue_index) for a legacy six-column table."""
    if requested == "legacy6-bitscore-evalue":
        return 4, 5
    if requested == "legacy6-evalue-bitscore":
        return 5, 4
    if requested not in {"auto", "legacy6-auto"}:
        raise SchemaError(f"{source}: unsupported six-column schema option {requested!r}")
    a = parse_float(values[4], "column 5", source)
    b = parse_float(values[5], "column 6", source)
    # Scientific notation / zero is strongly indicative of an E value.  A
    # numeric heuristic remains only for old unheadered files and deliberately
    # rejects ambiguous pairs rather than inventing a result.
    a_evalue_like = "e" in values[4].lower() or a == 0.0 or (0 < a < 1)
    b_evalue_like = "e" in values[5].lower() or b == 0.0 or (0 < b < 1)
    if a_evalue_like and not b_evalue_like:
        return 5, 4
    if b_evalue_like and not a_evalue_like:
        return 4, 5
    raise SchemaError(
        f"{source}: cannot safely infer six-column order from {values[4]!r}, {values[5]!r}; "
        "pass --blast-schema legacy6-bitscore-evalue or legacy6-evalue-bitscore."
    )


def iter_blast_hits(path: str | Path, schema: str = "auto") -> tuple[Iterator[BlastHit], str]:
    """Return an iterator of parsed hits and its detected schema label.

    Accepted inputs are standard unheadered 12-column outfmt 6, headered TSV
    with explicit names, or legacy six-column compact tables.  A ``#`` comment
    before data is ignored; a commented outfmt specification is not used as a
    schema because it is often stale in archived jobs.
    """
    source = str(path)
    with open(path, encoding="utf-8") as handle:
        first_values: list[str] | None = None
        first_line_number = 0
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            first_values = _split(raw)
            first_line_number = line_number
            break
    if first_values is None:
        raise SchemaError(f"{path}: no BLAST data rows")

    header_mapping: dict[str, int] | None = None
    if _looks_like_header(first_values):
        header_mapping = _canonical_header(first_values, f"{path}:{first_line_number}")
        detected = "headered"
    elif len(first_values) == 12:
        if schema not in {"auto", "blast12"}:
            raise SchemaError(f"{path}: detected 12 columns but --blast-schema={schema!r}")
        detected = "blast12"
    elif len(first_values) == 6:
        if schema not in {"auto", "legacy6-auto", "legacy6-bitscore-evalue", "legacy6-evalue-bitscore"}:
            raise SchemaError(f"{path}: detected 6 columns but --blast-schema={schema!r}")
        bitscore_index, evalue_index = _six_column_mapping(first_values, f"{path}:{first_line_number}", schema)
        detected = "legacy6-bitscore-evalue" if bitscore_index == 4 else "legacy6-evalue-bitscore"
    else:
        raise SchemaError(
            f"{path}:{first_line_number}: expected 12 standard BLAST columns or a documented 6-column legacy schema; "
            f"found {len(first_values)} columns"
        )

    def generator() -> Iterator[BlastHit]:
        with open(path, encoding="utf-8") as handle:
            skipped_header = header_mapping is not None
            for line_number, raw in enumerate(handle, start=1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                values = _split(raw)
                if skipped_header:
                    skipped_header = False
                    continue
                row_source = f"{path}:{line_number}"
                if header_mapping is not None:
                    if len(values) < len(header_mapping):
                        raise SchemaError(f"{row_source}: fewer fields than its header")
                    get = lambda name: values[header_mapping[name]]
                    def optional(name: str) -> int | None:
                        return parse_int(get(name), name, row_source) if name in header_mapping and get(name) else None
                    yield BlastHit(
                        get("qseqid"), get("sseqid"),
                        parse_float(get("pident"), "pident", row_source),
                        parse_int(get("length"), "length", row_source),
                        parse_float(get("evalue"), "evalue", row_source),
                        parse_float(get("bitscore"), "bitscore", row_source),
                        optional("qstart"), optional("qend"), optional("sstart"), optional("send"), line_number,
                    )
                elif detected == "blast12":
                    if len(values) != 12:
                        raise SchemaError(f"{row_source}: mixed BLAST schemas; expected 12 columns, found {len(values)}")
                    yield BlastHit(
                        values[0], values[1], parse_float(values[2], "pident", row_source),
                        parse_int(values[3], "length", row_source), parse_float(values[10], "evalue", row_source),
                        parse_float(values[11], "bitscore", row_source), parse_int(values[6], "qstart", row_source),
                        parse_int(values[7], "qend", row_source), parse_int(values[8], "sstart", row_source),
                        parse_int(values[9], "send", row_source), line_number,
                    )
                else:
                    if len(values) != 6:
                        raise SchemaError(f"{row_source}: mixed BLAST schemas; expected 6 columns, found {len(values)}")
                    yield BlastHit(
                        values[0], values[1], parse_float(values[2], "pident", row_source),
                        parse_int(values[3], "length", row_source), parse_float(values[evalue_index], "evalue", row_source),
                        parse_float(values[bitscore_index], "bitscore", row_source), raw_line_number=line_number,
                    )
    return generator(), detected


def _interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return max(start_a, start_b) <= min(end_a, end_b)


def _is_qualified(hit: BlastHit, min_identity: float, min_bitscore: float, max_evalue: float, min_alignment_length: int) -> bool:
    return (
        hit.pident >= min_identity and hit.bitscore >= min_bitscore and hit.length >= min_alignment_length
        and hit.evalue < max_evalue
    )


def _best_key(hit: BlastHit) -> tuple[float, float, float, int, str, int]:
    # `math.inf` handles the theoretical but invalid case of missing values.
    return (hit.evalue, -hit.bitscore, -hit.pident, -hit.length, hit.sseqid, hit.raw_line_number)


def classify_tblastx(
    candidate_path: str | Path,
    blast_path: str | Path,
    output_path: str | Path,
    schema_path: str | Path,
    blast_schema: str = "auto",
    min_identity: float = 50.0,
    min_bitscore: float = 50.0,
    max_evalue: float = 1e-5,
    min_alignment_length: int = 0,
    strictness: str = "legacy",
    synteny_padding_bp: int = 0,
    uncertain_ids_path: str | Path | None = None,
    query_fasta_path: str | Path | None = None,
    compatibility_lists_dir: str | Path | None = None,
) -> list[dict[str, object]]:
    """Classify candidate genes as pseudogenized, deleted, or uncertain.

    ``legacy`` reproduces the manuscript-era decision rule: any qualifying
    genome-wide tBLASTX hit means ``pseudogenized``.  ``synteny-aware`` requires
    candidate columns ``expected_target_chromosome``, ``expected_target_start``
    and ``expected_target_end`` and accepts only hits overlapping that interval
    (with optional padding).  It is the safer mode for revisions because it
    reduces remote-paralog/repeat false positives.
    """
    if strictness not in {"legacy", "synteny-aware"}:
        raise SchemaError("strictness must be 'legacy' or 'synteny-aware'")
    candidates = read_tsv(candidate_path, required=["target_sample", "reference_gene"])
    candidate_ids = [row["reference_gene"] for row in candidates]
    if len(set(candidate_ids)) != len(candidate_ids):
        raise SchemaError(f"{candidate_path}: reference_gene must be unique for a single target sample")
    target_samples = {row["target_sample"] for row in candidates}
    if len(target_samples) != 1:
        raise SchemaError(f"{candidate_path}: one target_sample per classification run is required")
    if strictness == "synteny-aware":
        needed = {"expected_target_chromosome", "expected_target_start", "expected_target_end"}
        missing = needed - set(candidates[0]) if candidates else needed
        if missing:
            raise SchemaError(
                f"{candidate_path}: synteny-aware mode requires columns {', '.join(sorted(missing))}; "
                "do not use a global tBLASTX hit as synteny evidence."
            )
    uncertain_ids = read_id_file(uncertain_ids_path) if uncertain_ids_path else set()
    query_fasta_ids = {identifier for identifier, _ in iter_fasta(query_fasta_path)} if query_fasta_path else None
    hit_iterator, detected_schema = iter_blast_hits(blast_path, schema=blast_schema)
    raw_by_query: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in hit_iterator:
        if hit.qseqid in set(candidate_ids):
            raw_by_query[hit.qseqid].append(hit)
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        query_id = candidate["reference_gene"]
        raw_hits = raw_by_query.get(query_id, [])
        qualifying = [
            hit for hit in raw_hits
            if _is_qualified(hit, min_identity, min_bitscore, max_evalue, min_alignment_length)
        ]
        syntenic: list[BlastHit] = []
        if strictness == "synteny-aware":
            expected_chrom = candidate["expected_target_chromosome"]
            expected_start = parse_int(candidate["expected_target_start"], "expected_target_start", str(candidate_path)) - synteny_padding_bp
            expected_end = parse_int(candidate["expected_target_end"], "expected_target_end", str(candidate_path)) + synteny_padding_bp
            syntenic = [
                hit for hit in qualifying
                if hit.subject_start is not None and hit.sseqid == expected_chrom
                and _interval_overlap(hit.subject_start, hit.subject_end or hit.subject_start, expected_start, expected_end)
            ]
            accepted = syntenic
        else:
            accepted = qualifying
        best = sorted(accepted, key=_best_key)[0] if accepted else None
        if query_fasta_ids is not None and query_id not in query_fasta_ids:
            # This is exactly what happened for eleven A. arguta A legacy
            # candidates: `sequences_extraction.py` removed IDs absent from
            # its reference nucleotide FASTA before BLAST.  A missing query is
            # missing evidence, not proof of a genomic deletion.
            classification, reason = "uncertain", "candidate_absent_from_tblastx_query_fasta"
        elif query_id in uncertain_ids:
            classification, reason = "uncertain", "listed_in_missing_data_or_gap_QC"
        elif best is not None:
            classification = "pseudogenized"
            reason = "qualifying_tBLASTX_hit" if strictness == "legacy" else "qualifying_hit_in_expected_syntenic_interval"
        elif strictness == "synteny-aware" and qualifying:
            classification, reason = "uncertain", "qualifying_global_hit_outside_expected_syntenic_interval"
        else:
            classification, reason = "deleted", "no_qualifying_tBLASTX_hit"
        rows.append({
            **candidate,
            "classification": classification,
            "decision_reason": reason,
            "tblastx_schema": detected_schema,
            "strictness": strictness,
            "n_raw_hits": len(raw_hits),
            "n_threshold_hits": len(qualifying),
            "n_syntenic_hits": len(syntenic) if strictness == "synteny-aware" else "",
            "best_subject_id": best.sseqid if best else "",
            "best_subject_start": best.subject_start if best and best.subject_start is not None else "",
            "best_subject_end": best.subject_end if best and best.subject_end is not None else "",
            "best_query_start": best.qstart if best and best.qstart is not None else "",
            "best_query_end": best.qend if best and best.qend is not None else "",
            "best_pident": format_number(best.pident) if best else "",
            "best_alignment_length": best.length if best else "",
            "best_bitscore": format_number(best.bitscore) if best else "",
            "best_evalue": format_number(best.evalue) if best else "",
            "best_raw_line": best.raw_line_number if best else "",
        })
    fields = list(candidates[0].keys()) if candidates else ["target_sample", "reference_gene"]
    fields += [
        "classification", "decision_reason", "tblastx_schema", "strictness", "n_raw_hits", "n_threshold_hits", "n_syntenic_hits",
        "best_subject_id", "best_subject_start", "best_subject_end", "best_query_start", "best_query_end",
        "best_pident", "best_alignment_length", "best_bitscore", "best_evalue", "best_raw_line",
    ]
    write_tsv(output_path, rows, fields)
    schema_rows = [{
        "blast_file": str(blast_path), "detected_schema": detected_schema, "min_identity": min_identity,
        "min_bitscore": min_bitscore, "max_evalue_exclusive": max_evalue,
        "min_alignment_length": min_alignment_length, "strictness": strictness,
        "candidate_file": str(candidate_path), "candidate_count": len(candidates),
        "query_fasta": str(query_fasta_path) if query_fasta_path else "",
        "candidate_ids_absent_from_query_fasta": sum(query_fasta_ids is not None and row["reference_gene"] not in query_fasta_ids for row in candidates),
    }]
    write_tsv(schema_path, schema_rows, [
        "blast_file", "detected_schema", "min_identity", "min_bitscore", "max_evalue_exclusive",
        "min_alignment_length", "strictness", "candidate_file", "candidate_count",
        "query_fasta", "candidate_ids_absent_from_query_fasta",
    ])
    if compatibility_lists_dir:
        output_dir = Path(compatibility_lists_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        sample = next(iter(target_samples))
        aliases = {"pseudogenized": "decayed", "deleted": "deleted", "uncertain": "uncertain"}
        for category, legacy_name in aliases.items():
            with open(output_dir / f"{sample}_{legacy_name}_genes.txt", "w", encoding="utf-8") as handle:
                for row in rows:
                    if row["classification"] == category:
                        handle.write(f"{row['reference_gene']}\n")
    return rows


def summarize_classification(
    classification_path: str | Path,
    reference_coords_path: str | Path,
    output_path: str | Path,
) -> list[dict[str, object]]:
    """Make an explicit per-sample loss summary from a classification table."""
    from .synorth import read_reference_coords

    rows = read_tsv(classification_path, required=["target_sample", "reference_gene", "classification"])
    sample_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sample_groups[row["target_sample"]].append(row)
    reference_gene_count = len(read_reference_coords(reference_coords_path))
    summary: list[dict[str, object]] = []
    for sample, group in sorted(sample_groups.items()):
        counts = {name: sum(row["classification"] == name for row in group) for name in ("pseudogenized", "deleted", "uncertain")}
        assessed = counts["pseudogenized"] + counts["deleted"]
        summary.append({
            "sample_id": sample,
            "reference_gene_count": reference_gene_count,
            "putative_candidate_count": len(group),
            "pseudogenized_count": counts["pseudogenized"],
            "deleted_count": counts["deleted"],
            "uncertain_count": counts["uncertain"],
            "assessed_loss_count": assessed,
            "putative_candidate_rate": len(group) / reference_gene_count,
            "assessed_loss_rate": assessed / reference_gene_count,
            "pseudogenized_rate": counts["pseudogenized"] / reference_gene_count,
            "deleted_rate": counts["deleted"] / reference_gene_count,
            "assessment_completion_rate": assessed / len(group) if group else 0.0,
        })
    write_tsv(output_path, summary, [
        "sample_id", "reference_gene_count", "putative_candidate_count", "pseudogenized_count", "deleted_count", "uncertain_count",
        "assessed_loss_count", "putative_candidate_rate", "assessed_loss_rate", "pseudogenized_rate", "deleted_rate", "assessment_completion_rate",
    ])
    return summary
