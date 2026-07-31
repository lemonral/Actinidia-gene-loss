"""Normalize SynOrths evidence and call *putative* loss candidates.

SynOrths itself remains an external program.  This module deliberately does
not hide its executable or parameters: the exact command is emitted by the
caller and the raw output is retained.  Candidate calls are not final losses;
they must be validated with the tBLASTX stage.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from .io_utils import SchemaError, natural_key, parse_float, parse_int, read_tsv, write_tsv


SYNORTH_FIELDS = [
    "first_gene", "first_chromosome", "first_start", "first_end",
    "second_gene", "second_chromosome", "second_start", "second_end",
    "evalue", "first_flank", "second_flank", "strand", "synorth_evidence",
]


def read_reference_coords(path: str | Path) -> list[dict[str, object]]:
    """Read either canonical headered coords TSV or legacy five-column coords."""
    source = str(path)
    with open(path, encoding="utf-8", newline="") as handle:
        first = handle.readline()
        if not first:
            raise SchemaError(f"{path}: empty coordinate file")
        handle.seek(0)
        header = first.rstrip("\n").split("\t")
        if "transcript_id" in header or "gene" in header:
            reader = csv.DictReader(handle, delimiter="\t")
            rows: list[dict[str, object]] = []
            for index, row in enumerate(reader, start=2):
                identifier = (row.get("transcript_id") or row.get("gene") or "").strip()
                chromosome = (row.get("chromosome") or row.get("chr") or "").strip()
                if not identifier or not chromosome:
                    raise SchemaError(f"{path}:{index}: coordinate row needs transcript_id/gene and chromosome/chr")
                rows.append({
                    "reference_gene": identifier,
                    "reference_chromosome": chromosome,
                    "reference_start": parse_int(row.get("start", ""), "start", f"{path}:{index}"),
                    "reference_end": parse_int(row.get("end", ""), "end", f"{path}:{index}"),
                    "reference_strand": (row.get("strand") or ".").strip(),
                })
        else:
            rows = []
            for index, raw in enumerate(handle, start=1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                values = raw.split()
                if len(values) < 5:
                    raise SchemaError(f"{path}:{index}: legacy coords needs >=5 whitespace-separated fields")
                rows.append({
                    "reference_gene": values[0], "reference_chromosome": values[1],
                    "reference_start": parse_int(values[2], "start", f"{path}:{index}"),
                    "reference_end": parse_int(values[3], "end", f"{path}:{index}"),
                    "reference_strand": values[4],
                })
    if not rows:
        raise SchemaError(f"{path}: no coordinate records")
    seen = set()
    for row in rows:
        identifier = str(row["reference_gene"])
        if identifier in seen:
            raise SchemaError(f"{path}: duplicate reference identifier {identifier!r}")
        seen.add(identifier)
    return sorted(rows, key=lambda row: (natural_key(str(row["reference_chromosome"])), int(row["reference_start"]), str(row["reference_gene"])))


def _read_raw_synorth(path: str | Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            values = line.split()
            if len(values) < 13:
                raise SchemaError(
                    f"{path}:{line_number}: SynOrths row needs at least 13 columns; found {len(values)}"
                )
            records.append({field: values[index] for index, field in enumerate(SYNORTH_FIELDS)})
    if not records:
        raise SchemaError(f"{path}: no SynOrths records")
    return records


def _detect_reference_side(records: list[dict[str, str]], reference_ids: set[str], requested: str) -> str:
    if requested in {"first", "second"}:
        return requested
    first_matches = sum(record["first_gene"] in reference_ids for record in records)
    second_matches = sum(record["second_gene"] in reference_ids for record in records)
    if first_matches == second_matches:
        raise SchemaError(
            "Cannot automatically identify the SynOrths reference side: "
            f"first matches={first_matches}, second matches={second_matches}. "
            "Pass --reference-side first or --reference-side second after checking IDs."
        )
    return "first" if first_matches > second_matches else "second"


def normalize_synorth(
    raw_path: str | Path,
    reference_coords_path: str | Path,
    target_sample: str,
    output_path: str | Path,
    reference_side: str = "auto",
) -> tuple[list[dict[str, object]], str]:
    """Convert a raw SynOrths table into a documented, stable schema."""
    reference_coords = read_reference_coords(reference_coords_path)
    reference_ids = {str(row["reference_gene"]) for row in reference_coords}
    raw_records = _read_raw_synorth(raw_path)
    side = _detect_reference_side(raw_records, reference_ids, reference_side)
    standardized: list[dict[str, object]] = []
    unmatched_reference_ids = 0
    for record in raw_records:
        if side == "first":
            prefix, other = "first", "second"
        else:
            prefix, other = "second", "first"
        reference_gene = record[f"{prefix}_gene"]
        if reference_gene not in reference_ids:
            unmatched_reference_ids += 1
        standardized.append({
            "target_sample": target_sample,
            "reference_gene": reference_gene,
            "reference_chromosome": record[f"{prefix}_chromosome"],
            "reference_start": parse_int(record[f"{prefix}_start"], "reference_start", str(raw_path)),
            "reference_end": parse_int(record[f"{prefix}_end"], "reference_end", str(raw_path)),
            "target_gene": record[f"{other}_gene"],
            "target_chromosome": record[f"{other}_chromosome"],
            "target_start": parse_int(record[f"{other}_start"], "target_start", str(raw_path)),
            "target_end": parse_int(record[f"{other}_end"], "target_end", str(raw_path)),
            "evalue": record["evalue"],
            "reference_flank": record[f"{prefix}_flank"],
            "target_flank": record[f"{other}_flank"],
            "strand": record["strand"],
            "synorth_evidence": record["synorth_evidence"],
            "raw_synorth_file": str(raw_path),
        })
    if unmatched_reference_ids:
        raise SchemaError(
            f"{raw_path}: {unmatched_reference_ids}/{len(standardized)} selected reference IDs are absent "
            f"from {reference_coords_path}; verify --reference-side and identifier normalization."
        )
    standardized.sort(key=lambda row: (
        natural_key(str(row["reference_chromosome"])), int(row["reference_start"]), str(row["reference_gene"]),
        parse_float(str(row["evalue"]), "evalue", str(raw_path)),
    ))
    write_tsv(output_path, standardized, [
        "target_sample", "reference_gene", "reference_chromosome", "reference_start", "reference_end",
        "target_gene", "target_chromosome", "target_start", "target_end", "evalue",
        "reference_flank", "target_flank", "strand", "synorth_evidence", "raw_synorth_file",
    ])
    return standardized, side


def read_normalized_synorth(path: str | Path) -> list[dict[str, str]]:
    return read_tsv(path, required=["target_sample", "reference_gene", "reference_chromosome"])


def call_candidates(
    reference_coords_path: str | Path,
    normalized_synorth_path: str | Path,
    output_path: str | Path,
    flank_genes: int = 20,
    mode: str = "bracketed",
    min_anchors_each_side: int = 1,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Call putative reference-gene losses from SynOrths anchor absence.

    ``legacy-neighbor`` exactly reproduces the archival
    ``lost_gene_identification_new.py`` set rule.  ``bracketed`` is the
    recommended default: an unmatched reference gene must be bounded by
    observed anchors within the configurable neighbourhood on *both* sides.
    Neither mode converts a candidate into a deletion; that happens only after
    remnant-search classification.
    """
    if flank_genes < 1:
        raise SchemaError("flank_genes must be at least 1")
    if min_anchors_each_side < 1:
        raise SchemaError("min_anchors_each_side must be at least 1")
    if mode not in {"legacy-neighbor", "bracketed", "all-unmatched"}:
        raise SchemaError("mode must be one of: legacy-neighbor, bracketed, all-unmatched")
    coords = read_reference_coords(reference_coords_path)
    synorth_rows = read_normalized_synorth(normalized_synorth_path)
    target_samples = {row["target_sample"] for row in synorth_rows}
    if len(target_samples) != 1:
        raise SchemaError(f"{normalized_synorth_path}: expected one target_sample, found {sorted(target_samples)}")
    target_sample = next(iter(target_samples))
    reference_ids = {str(row["reference_gene"]) for row in coords}
    anchors = {row["reference_gene"] for row in synorth_rows}
    unknown = anchors - reference_ids
    if unknown:
        preview = ", ".join(sorted(unknown)[:5])
        raise SchemaError(f"{normalized_synorth_path}: SynOrths contains unknown reference IDs (e.g. {preview})")

    by_chromosome: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in coords:
        by_chromosome[str(row["reference_chromosome"])].append(row)

    candidates: list[dict[str, object]] = []
    retained: list[dict[str, object]] = []
    for chromosome in sorted(by_chromosome, key=natural_key):
        records = by_chromosome[chromosome]
        ordered_ids = [str(row["reference_gene"]) for row in records]
        anchor_positions = {index for index, identifier in enumerate(ordered_ids) if identifier in anchors}
        for index, row in enumerate(records):
            identifier = str(row["reference_gene"])
            left_positions = range(max(0, index - flank_genes), index)
            right_positions = range(index + 1, min(len(records), index + flank_genes + 1))
            left_count = sum(position in anchor_positions for position in left_positions)
            right_count = sum(position in anchor_positions for position in right_positions)
            common = {
                "target_sample": target_sample,
                "reference_gene": identifier,
                "reference_chromosome": chromosome,
                "reference_start": row["reference_start"],
                "reference_end": row["reference_end"],
                "reference_strand": row["reference_strand"],
                "flank_genes": flank_genes,
                "left_anchor_count": left_count,
                "right_anchor_count": right_count,
            }
            if identifier in anchors:
                retained.append({**common, "synorth_status": "retained_anchor", "candidate_rule": ""})
                continue
            if mode == "legacy-neighbor":
                # Exact old logic: a gene is a candidate if it occurs in the
                # window around at least one anchor; anchors too close to either
                # chromosome end contribute no window at all.
                legacy_anchor_positions = {
                    position for position in anchor_positions
                    if position >= flank_genes and position <= len(records) - flank_genes - 1
                }
                is_candidate = any(abs(index - position) <= flank_genes and index != position for position in legacy_anchor_positions)
                rule = "legacy_neighbor_of_non_edge_anchor"
            elif mode == "bracketed":
                is_candidate = left_count >= min_anchors_each_side and right_count >= min_anchors_each_side
                rule = "unmatched_bracketed_by_anchor_windows"
            else:
                is_candidate = True
                rule = "all_unmatched_reference_genes"
            if is_candidate:
                candidates.append({**common, "synorth_status": "putative_loss", "candidate_rule": rule})
    candidate_fields = [
        "target_sample", "reference_gene", "reference_chromosome", "reference_start", "reference_end", "reference_strand",
        "synorth_status", "candidate_rule", "flank_genes", "left_anchor_count", "right_anchor_count",
    ]
    candidates.sort(key=lambda row: (natural_key(str(row["reference_chromosome"])), int(row["reference_start"]), str(row["reference_gene"])))
    retained.sort(key=lambda row: (natural_key(str(row["reference_chromosome"])), int(row["reference_start"]), str(row["reference_gene"])))
    write_tsv(output_path, candidates, candidate_fields)
    retained_path = Path(output_path).with_name(Path(output_path).stem + ".retained_anchors.tsv")
    write_tsv(retained_path, retained, candidate_fields)
    return candidates, retained
