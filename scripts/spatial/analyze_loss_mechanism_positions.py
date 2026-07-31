#!/usr/bin/env python3
"""Classify loss evidence and analyze residual-sequence chromosome positions.

The loss numerator is the frozen decayed-plus-deleted classification.  Existing
Miniprot output is used only as an orthogonal refinement: explicit disruptions,
partial/truncated alignments, unresolved residual sequence, and absence of the
historical translated-search signal remain mutually exclusive.  Local
alignments retain their observed target-genome coordinates.  For calls without
local support, the best already-generated genome-wide Miniprot alignment is
reported as a candidate placement; a different-chromosome alignment is never
described as a proven rearrangement because it can also represent a paralog.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


LOSS_TYPES = (
    "no_qualifying_translated_hit",
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
    "truncation_or_partial_alignment_candidate",
    "residual_sequence_mechanism_unresolved",
)
STRICT_TYPES = (
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
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
RELATIONS = (
    "expected_interval_local",
    "same_chromosome_displacement_candidate",
    "interchromosomal_displacement_candidate",
    "genomewide_residual_sequence_unanchored",
    "unlocalized",
)


class SpatialMechanismError(ValueError):
    """Raised when frozen evidence cannot be reconciled without inference."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article-matrix", required=True, type=Path)
    parser.add_argument("--uniform-config", required=True, type=Path)
    parser.add_argument("--chromosome-map-config", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--reference-gff", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-chromosomes", type=int, default=29)
    parser.add_argument("--expected-positive-rows", type=int, default=179827)
    parser.add_argument("--minimum-query-coverage", type=float, default=0.50)
    parser.add_argument("--minimum-identity", type=float, default=0.50)
    parser.add_argument("--minimum-alignment-score", type=int, default=50)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    if not path.is_file() or (path.stat().st_size == 0 and not allow_empty):
        raise SpatialMechanismError(f"missing or empty input: {path}")
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
        raise SpatialMechanismError(f"invalid TSV: {path.name}")
    return rows, fields


def require(fields: list[str], required: set[str], label: str) -> None:
    missing = sorted(required - set(fields))
    if missing:
        raise SpatialMechanismError(f"{label} missing fields: {', '.join(missing)}")


def resolve(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise SpatialMechanismError(f"unsafe data-root-relative path: {value!r}")
    result = (root / relative).absolute()
    if not result.is_relative_to(root):
        raise SpatialMechanismError(f"path escapes data root: {value!r}")
    return result


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


def loss_type(cause: str) -> str:
    if cause == "no_qualifying_genomewide_tblastx_hit":
        return "no_qualifying_translated_hit"
    if cause == "frameshift_supported":
        return "frameshift_supported"
    if cause == "stop_supported":
        return "inframe_stop_supported"
    if cause == "frameshift_and_stop_supported":
        return "frameshift_and_stop_supported"
    if cause in TRUNCATION_CAUSES:
        return "truncation_or_partial_alignment_candidate"
    if cause in UNRESOLVED_CAUSES:
        return "residual_sequence_mechanism_unresolved"
    raise SpatialMechanismError(f"unmapped positive-loss cause: {cause!r}")


def integer(value: str, context: str) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise SpatialMechanismError(f"{context}: invalid integer {value!r}") from error


def parse_paf(raw: str, context: str) -> dict[str, object]:
    fields = raw.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise SpatialMechanismError(f"{context}: fewer than 12 PAF fields")
    qlen = integer(fields[1], context)
    qstart = integer(fields[2], context)
    qend = integer(fields[3], context)
    target_length = integer(fields[6], context)
    target_start0 = integer(fields[7], context)
    target_end0 = integer(fields[8], context)
    matches = integer(fields[9], context)
    aligned = integer(fields[10], context)
    if (
        qlen < 1
        or not 0 <= qstart < qend <= qlen
        or target_length < 1
        or not 0 <= target_start0 < target_end0 <= target_length
        or aligned < 1
        or not 0 <= matches <= aligned
        or fields[4] not in {"+", "-"}
    ):
        raise SpatialMechanismError(f"{context}: invalid PAF coordinates")
    tags: dict[str, str] = {}
    for item in fields[12:]:
        pieces = item.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]
    for tag in ("AS", "fs", "st"):
        if tag not in tags:
            raise SpatialMechanismError(f"{context}: missing PAF {tag} tag")
    return {
        "query": fields[0],
        "query_length": qlen,
        "query_start": qstart,
        "query_end": qend,
        "strand": fields[4],
        "target": fields[5],
        "target_length": target_length,
        "target_start": target_start0 + 1,
        "target_end": target_end0,
        "coverage": (qend - qstart) / qlen,
        "identity": matches / aligned,
        "score": integer(tags["AS"], context),
        "frameshifts": integer(tags["fs"], context),
        "stops": integer(tags["st"], context),
    }


def better_alignment(
    candidate: dict[str, object],
    current: dict[str, object] | None,
) -> bool:
    if current is None:
        return True
    candidate_key = (
        int(candidate["score"]),
        float(candidate["coverage"]),
        float(candidate["identity"]),
        -int(candidate["target_start"]),
    )
    current_key = (
        int(current["score"]),
        float(current["coverage"]),
        float(current["identity"]),
        -int(current["target_start"]),
    )
    return candidate_key > current_key


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def chi_square_survival(value: float, degrees_of_freedom: int) -> float:
    """Return chi-square survival probability without an optional dependency."""
    if value < 0 or degrees_of_freedom < 1:
        raise SpatialMechanismError("invalid chi-square statistic")
    a = degrees_of_freedom / 2.0
    x = value / 2.0
    if x == 0:
        return 1.0
    epsilon = 3e-14
    floor = 1e-300
    maximum_iterations = 10000
    if x < a + 1.0:
        # Regularized lower incomplete gamma by its convergent series.
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(maximum_iterations):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * epsilon:
                lower = total * math.exp(-x + a * math.log(x) - math.lgamma(a))
                return min(max(1.0 - lower, 0.0), 1.0)
        raise SpatialMechanismError("incomplete-gamma series did not converge")
    # Regularized upper incomplete gamma by Lentz's continued fraction.
    b = x + 1.0 - a
    c = 1.0 / floor
    d = 1.0 / max(b, floor)
    fraction = d
    for iteration in range(1, maximum_iterations + 1):
        coefficient = -iteration * (iteration - a)
        b += 2.0
        d = coefficient * d + b
        if abs(d) < floor:
            d = floor
        c = b + coefficient / c
        if abs(c) < floor:
            c = floor
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) < epsilon:
            upper = math.exp(-x + a * math.log(x) - math.lgamma(a)) * fraction
            return min(max(upper, 0.0), 1.0)
    raise SpatialMechanismError("incomplete-gamma fraction did not converge")


def bh_adjust(rows: list[dict[str, object]], key: str, output: str) -> None:
    ordered = sorted(
        enumerate(rows),
        key=lambda item: float(item[1][key]),
    )
    count = len(ordered)
    running = 1.0
    for reverse_rank, (index, row) in enumerate(reversed(ordered), 1):
        rank = count - reverse_rank + 1
        adjusted = min(running, float(row[key]) * count / rank, 1.0)
        running = adjusted
        rows[index][output] = f"{adjusted:.12g}"


def audit_te_annotation(path: Path) -> tuple[dict[str, object], Counter[str]]:
    feature_counts: Counter[str] = Counter()
    repeat_like = 0
    keywords = (
        "repeat",
        "transpos",
        "retro",
        "ltr",
        "line",
        "sine",
        "gypsy",
        "copia",
    )
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SpatialMechanismError(
                    f"{path.name}:{line_number}: GFF row does not have 9 fields"
                )
            feature_counts[fields[2]] += 1
            if any(keyword in raw.lower() for keyword in keywords):
                repeat_like += 1
    if not feature_counts:
        raise SpatialMechanismError("reference GFF contains no features")
    status = (
        "AVAILABLE_REPEAT_FEATURES_PRESENT"
        if repeat_like
        else "UNAVAILABLE_NO_TE_ANNOTATION"
    )
    return (
        {
            "reference_annotation": path.name,
            "status": status,
            "repeat_or_transposon_feature_rows": repeat_like,
            "te_association_performed": "false",
            "reason": (
                "TE association requires an independently supplied repeat annotation"
                if not repeat_like
                else "repeat-like rows require a dedicated validated TE adapter"
            ),
        },
        feature_counts,
    )


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise SpatialMechanismError(
            f"refusing to overwrite output directory: {args.output_dir}"
        )
    if not (
        0 < args.minimum_query_coverage <= 1
        and 0 < args.minimum_identity <= 1
        and args.minimum_alignment_score >= 0
    ):
        raise SpatialMechanismError("invalid alignment thresholds")
    root = args.data_root.resolve()
    inputs = [
        args.article_matrix.resolve(),
        args.uniform_config.resolve(),
        args.chromosome_map_config.resolve(),
        args.unit_metadata.resolve(),
        args.reference_gff.resolve(),
    ]
    for path in inputs:
        binding(path)

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
    metadata = {
        row["assembly_unit_id"]: row
        for row in metadata_rows
        if row["include"].lower() == "true"
    }
    if len(metadata) != args.expected_units:
        raise SpatialMechanismError("unit metadata count changed")

    uniform_rows, uniform_fields = read_tsv(args.uniform_config)
    require(
        uniform_fields,
        {"unit", "target_genome", "output_dir"},
        args.uniform_config.name,
    )
    uniform = {row["unit"]: row for row in uniform_rows}
    if set(uniform) != set(metadata) or len(uniform) != len(uniform_rows):
        raise SpatialMechanismError("uniform config does not match unit metadata")

    map_rows, map_fields = read_tsv(args.chromosome_map_config)
    require(
        map_fields,
        {"unit", "map_path", "map_mode"},
        args.chromosome_map_config.name,
    )
    map_config = {row["unit"]: row for row in map_rows}
    if set(map_config) != set(metadata) or len(map_config) != len(map_rows):
        raise SpatialMechanismError("chromosome-map config does not match unit metadata")

    positives: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    positive_rows = 0
    with open_text(args.article_matrix) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        require(
            fields,
            {
                "reference_gene_id",
                "assembly_unit_id",
                "source_group",
                "manuscript_classification",
                "manuscript_positive_loss",
                "callable",
                "refined_decayed_cause",
                "refined_cause_evidence_level",
            },
            args.article_matrix.name,
        )
        for line_number, row in enumerate(reader, 2):
            if row["manuscript_positive_loss"] != "true":
                continue
            unit = row["assembly_unit_id"]
            gene = row["reference_gene_id"]
            if unit not in metadata or gene in positives[unit]:
                raise SpatialMechanismError(
                    f"{args.article_matrix.name}:{line_number}: invalid positive key"
                )
            row["loss_type_group"] = loss_type(row["refined_decayed_cause"])
            positives[unit][gene] = row
            positive_rows += 1
    if args.expected_positive_rows and positive_rows != args.expected_positive_rows:
        raise SpatialMechanismError(
            f"positive loss row count changed: {positive_rows}"
        )

    detail_rows: list[dict[str, object]] = []
    chromosome_lengths: dict[tuple[str, str], int] = {}
    unit_input_paths: list[Path] = []
    for unit in sorted(metadata):
        mapping: dict[str, str] = {}
        mapping_row = map_config[unit]
        if mapping_row["map_mode"] == "already_hy4a_chr_labels":
            mapping = {
                f"Chr{index:02d}": f"Chr{index:02d}"
                for index in range(1, args.expected_chromosomes + 1)
            }
        elif mapping_row["map_mode"] == "similarity_map_to_hy4a":
            map_path = resolve(root, mapping_row["map_path"])
            rows, fields = read_tsv(map_path)
            require(
                fields,
                {"query_chromosome", "final_chromosome"},
                map_path.name,
            )
            mapping = {
                row["query_chromosome"]: row["final_chromosome"]
                for row in rows
            }
            unit_input_paths.append(map_path)
        else:
            raise SpatialMechanismError(
                f"{unit}: unknown map mode {mapping_row['map_mode']!r}"
            )
        if len(set(mapping.values())) != args.expected_chromosomes:
            raise SpatialMechanismError(
                f"{unit}: chromosome map is not a {args.expected_chromosomes}-chromosome set"
            )

        output_dir = resolve(root, uniform[unit]["output_dir"])
        state_path = output_dir / "uniform_candidate_loss_states.tsv"
        paf_path = output_dir / "raw_alignments.paf.gz"
        state_rows, state_fields = read_tsv(state_path)
        require(
            state_fields,
            {
                "reference_gene",
                "callable",
                "target_chromosome",
                "target_interval_start_1based",
                "target_interval_end_1based",
                "qualifying_local_alignment",
                "alignment_target_start_1based",
                "alignment_target_end_1based",
                "alignment_strand",
                "query_coverage",
                "exact_alignment_identity",
                "alignment_score",
                "frameshift_events",
                "inframe_stop_codons",
            },
            state_path.name,
        )
        state = {
            row["reference_gene"]: row
            for row in state_rows
            if row["reference_gene"] in positives[unit]
        }
        unit_input_paths.extend([state_path, paf_path])

        best: dict[str, dict[str, object]] = {}
        with gzip.open(paf_path, "rt", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                alignment = parse_paf(raw, f"{unit}:PAF:{line_number}")
                raw_chromosome = str(alignment["target"])
                if raw_chromosome in mapping:
                    final_chromosome = mapping[raw_chromosome]
                    length_key = (unit, final_chromosome)
                    length = int(alignment["target_length"])
                    if (
                        length_key in chromosome_lengths
                        and chromosome_lengths[length_key] != length
                    ):
                        raise SpatialMechanismError(
                            f"{unit}: inconsistent length for {final_chromosome}"
                        )
                    chromosome_lengths[length_key] = length
                query = str(alignment["query"])
                if query not in positives[unit]:
                    continue
                qualifying = (
                    float(alignment["coverage"]) >= args.minimum_query_coverage
                    and float(alignment["identity"]) >= args.minimum_identity
                    and int(alignment["score"]) >= args.minimum_alignment_score
                )
                if qualifying and better_alignment(alignment, best.get(query)):
                    best[query] = alignment

        expected_chromosome_set = {
            f"Chr{index:02d}"
            for index in range(1, args.expected_chromosomes + 1)
        }
        observed_lengths = {
            chromosome
            for observed_unit, chromosome in chromosome_lengths
            if observed_unit == unit
        }
        if observed_lengths != expected_chromosome_set:
            missing = sorted(expected_chromosome_set - observed_lengths)
            raise SpatialMechanismError(
                f"{unit}: PAF lacks target lengths for {', '.join(missing)}"
            )

        for gene in sorted(positives[unit]):
            evidence = positives[unit][gene]
            state_row = state.get(gene)
            expected_raw = ""
            expected_final = ""
            expected_start = ""
            expected_end = ""
            local = False
            selected: dict[str, object] | None = None
            position_source = ""
            if state_row is not None:
                expected_raw = state_row["target_chromosome"]
                expected_final = mapping.get(expected_raw, "")
                expected_start = state_row["target_interval_start_1based"]
                expected_end = state_row["target_interval_end_1based"]
                local = state_row["qualifying_local_alignment"] == "true"
                if local:
                    selected = {
                        "target": expected_raw,
                        "target_length": chromosome_lengths[(unit, expected_final)],
                        "target_start": integer(
                            state_row["alignment_target_start_1based"],
                            f"{unit}:{gene}:local start",
                        ),
                        "target_end": integer(
                            state_row["alignment_target_end_1based"],
                            f"{unit}:{gene}:local end",
                        ),
                        "strand": state_row["alignment_strand"],
                        "coverage": float(state_row["query_coverage"]),
                        "identity": float(state_row["exact_alignment_identity"]),
                        "score": integer(
                            state_row["alignment_score"],
                            f"{unit}:{gene}:local score",
                        ),
                        "frameshifts": integer(
                            state_row["frameshift_events"],
                            f"{unit}:{gene}:local frameshift",
                        ),
                        "stops": integer(
                            state_row["inframe_stop_codons"],
                            f"{unit}:{gene}:local stop",
                        ),
                    }
                    position_source = "qualifying_local_miniprot_alignment"
            if selected is None and gene in best:
                selected = best[gene]
                position_source = "best_qualifying_genomewide_miniprot_alignment"

            relation = "unlocalized"
            residual_raw = ""
            residual_final = ""
            residual_start: int | str = ""
            residual_end: int | str = ""
            residual_midpoint: float | str = ""
            normalized_end_distance: float | str = ""
            target_length: int | str = ""
            strand = ""
            coverage: float | str = ""
            identity: float | str = ""
            score: int | str = ""
            frameshifts: int | str = ""
            stops: int | str = ""
            if selected is not None:
                residual_raw = str(selected["target"])
                residual_final = mapping.get(residual_raw, "")
                if not residual_final:
                    relation = "unlocalized"
                    position_source = (
                        "qualifying_alignment_on_nonprimary_sequence"
                    )
                else:
                    residual_start = int(selected["target_start"])
                    residual_end = int(selected["target_end"])
                    residual_midpoint = (residual_start + residual_end) / 2.0
                    target_length = chromosome_lengths[(unit, residual_final)]
                    normalized_end_distance = min(
                        residual_midpoint - 1.0,
                        target_length - residual_midpoint,
                    ) / ((target_length - 1.0) / 2.0)
                    strand = str(selected["strand"])
                    coverage = float(selected["coverage"])
                    identity = float(selected["identity"])
                    score = int(selected["score"])
                    frameshifts = int(selected["frameshifts"])
                    stops = int(selected["stops"])
                    if local:
                        relation = "expected_interval_local"
                    elif not expected_final:
                        relation = "genomewide_residual_sequence_unanchored"
                    elif residual_final == expected_final:
                        relation = "same_chromosome_displacement_candidate"
                    else:
                        relation = "interchromosomal_displacement_candidate"

            detail_rows.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": metadata[unit]["biological_species"],
                    "haplotype_or_subgenome": metadata[unit][
                        "haplotype_or_subgenome"
                    ],
                    "reference_gene_id": gene,
                    "source_group": evidence["source_group"],
                    "primary_classification": evidence[
                        "manuscript_classification"
                    ],
                    "refined_cause": evidence["refined_decayed_cause"],
                    "refined_cause_evidence_level": evidence[
                        "refined_cause_evidence_level"
                    ],
                    "loss_type_group": evidence["loss_type_group"],
                    "callable": evidence["callable"],
                    "expected_chromosome_raw": expected_raw,
                    "expected_chromosome_hy4a": expected_final,
                    "expected_interval_start_1based": expected_start,
                    "expected_interval_end_1based": expected_end,
                    "position_source": position_source,
                    "residual_chromosome_raw": residual_raw,
                    "residual_chromosome_hy4a": residual_final,
                    "residual_start_1based": residual_start,
                    "residual_end_1based": residual_end,
                    "residual_midpoint_1based": residual_midpoint,
                    "residual_strand": strand,
                    "target_chromosome_length": target_length,
                    "normalized_end_distance": normalized_end_distance,
                    "alignment_query_coverage": coverage,
                    "alignment_exact_identity": identity,
                    "alignment_score": score,
                    "alignment_frameshifts": frameshifts,
                    "alignment_inframe_stops": stops,
                    "location_relation": relation,
                    "spatial_eligible": str(relation != "unlocalized").lower(),
                }
            )

    if len(detail_rows) != positive_rows:
        raise SpatialMechanismError("detailed positive rows do not close")

    mechanism_rows: list[dict[str, object]] = []
    relation_rows: list[dict[str, object]] = []
    unit_mechanism_rows: list[dict[str, object]] = []
    for grouped_type in LOSS_TYPES:
        selected = [
            row for row in detail_rows if row["loss_type_group"] == grouped_type
        ]
        relation_counts = Counter(str(row["location_relation"]) for row in selected)
        distances = [
            float(row["normalized_end_distance"])
            for row in selected
            if row["normalized_end_distance"] != ""
        ]
        mechanism_rows.append(
            {
                "loss_type_group": grouped_type,
                "positive_unit_gene_rows": len(selected),
                "spatially_placed_rows": len(distances),
                "unlocalized_rows": relation_counts["unlocalized"],
                "median_normalized_end_distance": (
                    f"{quantile(distances, 0.5):.12g}" if distances else ""
                ),
                "q1_normalized_end_distance": (
                    f"{quantile(distances, 0.25):.12g}" if distances else ""
                ),
                "q3_normalized_end_distance": (
                    f"{quantile(distances, 0.75):.12g}" if distances else ""
                ),
            }
        )
        for relation in RELATIONS:
            relation_rows.append(
                {
                    "loss_type_group": grouped_type,
                    "location_relation": relation,
                    "unit_gene_rows": relation_counts[relation],
                }
            )
    for unit in sorted(metadata):
        selected = [row for row in detail_rows if row["assembly_unit_id"] == unit]
        counts = Counter(str(row["loss_type_group"]) for row in selected)
        for grouped_type in LOSS_TYPES:
            unit_mechanism_rows.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": metadata[unit]["biological_species"],
                    "haplotype_or_subgenome": metadata[unit][
                        "haplotype_or_subgenome"
                    ],
                    "loss_type_group": grouped_type,
                    "positive_unit_gene_rows": counts[grouped_type],
                }
            )

    chromosomes = [
        f"Chr{index:02d}" for index in range(1, args.expected_chromosomes + 1)
    ]
    chromosome_rows: list[dict[str, object]] = []
    heterogeneity_rows: list[dict[str, object]] = []
    for grouped_type in LOSS_TYPES:
        type_rows = [
            row
            for row in detail_rows
            if row["loss_type_group"] == grouped_type
            and row["spatial_eligible"] == "true"
        ]
        observed = Counter(str(row["residual_chromosome_hy4a"]) for row in type_rows)
        expected: Counter[str] = Counter()
        unit_totals = Counter(str(row["assembly_unit_id"]) for row in type_rows)
        for unit, total in unit_totals.items():
            unit_length_total = sum(
                chromosome_lengths[(unit, chromosome)]
                for chromosome in chromosomes
            )
            for chromosome in chromosomes:
                expected[chromosome] += (
                    total
                    * chromosome_lengths[(unit, chromosome)]
                    / unit_length_total
                )
        chi_square = 0.0
        for chromosome in chromosomes:
            if expected[chromosome] > 0:
                chi_square += (
                    observed[chromosome] - expected[chromosome]
                ) ** 2 / expected[chromosome]
        p_value = (
            chi_square_survival(chi_square, len(chromosomes) - 1)
            if type_rows
            else 1.0
        )
        heterogeneity_rows.append(
            {
                "loss_type_group": grouped_type,
                "spatially_placed_rows": len(type_rows),
                "chi_square": f"{chi_square:.12g}",
                "degrees_of_freedom": len(chromosomes) - 1,
                "p_value": f"{p_value:.12g}",
                "test_definition": (
                    "observed HY4A-standardized chromosome counts versus "
                    "within-unit chromosome-length opportunities"
                ),
            }
        )
        for chromosome in chromosomes:
            total_length = sum(
                chromosome_lengths[(unit, chromosome)]
                for unit in metadata
            )
            count = observed[chromosome]
            chromosome_rows.append(
                {
                    "loss_type_group": grouped_type,
                    "chromosome_hy4a": chromosome,
                    "observed_residual_rows": count,
                    "summed_chromosome_length_bp_across_units": total_length,
                    "residual_rows_per_100mb": f"{count * 1e8 / total_length:.12g}",
                    "length_opportunity_expected_rows": f"{expected[chromosome]:.12g}",
                }
            )
    bh_adjust(heterogeneity_rows, "p_value", "bh_q_value")

    strict_rows: list[dict[str, object]] = []
    for grouped_type in STRICT_TYPES:
        values = [
            float(row["normalized_end_distance"])
            for row in detail_rows
            if row["loss_type_group"] == grouped_type
            and row["normalized_end_distance"] != ""
        ]
        strict_rows.append(
            {
                "loss_type_group": grouped_type,
                "spatially_placed_rows": len(values),
                "minimum": f"{min(values):.12g}" if values else "",
                "q1": f"{quantile(values, 0.25):.12g}" if values else "",
                "median": f"{quantile(values, 0.5):.12g}" if values else "",
                "q3": f"{quantile(values, 0.75):.12g}" if values else "",
                "maximum": f"{max(values):.12g}" if values else "",
            }
        )

    te_row, feature_counts = audit_te_annotation(args.reference_gff)
    te_row["feature_type_counts"] = ";".join(
        f"{name}:{count}" for name, count in sorted(feature_counts.items())
    )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        detail_path = staging / "loss_residual_positions.tsv.gz"
        mechanism_path = staging / "loss_mechanism_summary.tsv"
        relation_path = staging / "location_relation_summary.tsv"
        unit_path = staging / "unit_loss_mechanism_summary.tsv"
        chromosome_path = staging / "chromosome_loss_mechanism_summary.tsv"
        heterogeneity_path = staging / "chromosome_heterogeneity_tests.tsv"
        strict_path = staging / "strict_loss_end_distance_summary.tsv"
        te_path = staging / "reference_te_annotation_audit.tsv"
        write_tsv_gz(detail_path, list(detail_rows[0]), detail_rows)
        write_tsv(mechanism_path, list(mechanism_rows[0]), mechanism_rows)
        write_tsv(relation_path, list(relation_rows[0]), relation_rows)
        write_tsv(unit_path, list(unit_mechanism_rows[0]), unit_mechanism_rows)
        write_tsv(chromosome_path, list(chromosome_rows[0]), chromosome_rows)
        write_tsv(
            heterogeneity_path,
            list(heterogeneity_rows[0]),
            heterogeneity_rows,
        )
        write_tsv(strict_path, list(strict_rows[0]), strict_rows)
        write_tsv(te_path, list(te_row), [te_row])
        output_paths = [
            detail_path,
            mechanism_path,
            relation_path,
            unit_path,
            chromosome_path,
            heterogeneity_path,
            strict_path,
            te_path,
        ]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_LOSS_MECHANISM_SPATIAL_ANALYSIS",
            "assembly_units": len(metadata),
            "hy4a_standardized_chromosomes": len(chromosomes),
            "positive_unit_gene_rows": positive_rows,
            "spatially_placed_rows": sum(
                row["spatial_eligible"] == "true" for row in detail_rows
            ),
            "unlocalized_rows": sum(
                row["spatial_eligible"] == "false" for row in detail_rows
            ),
            "loss_type_groups": list(LOSS_TYPES),
            "strict_disruption_types": list(STRICT_TYPES),
            "parameters": {
                "minimum_query_coverage": args.minimum_query_coverage,
                "minimum_exact_alignment_identity": args.minimum_identity,
                "minimum_alignment_score": args.minimum_alignment_score,
                "candidate_displacement_interpretation": (
                    "best existing genome-wide Miniprot alignment; not proof of "
                    "inversion, translocation, or orthology"
                ),
            },
            "te_annotation_status": te_row["status"],
            "inputs": [
                binding(path)
                for path in inputs + sorted(set(unit_input_paths))
            ],
            "outputs": [binding(path, allow_empty=True) for path in output_paths],
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
    except (
        SpatialMechanismError,
        OSError,
        UnicodeError,
        csv.Error,
        ValueError,
    ) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
