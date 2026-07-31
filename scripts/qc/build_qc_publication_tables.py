#!/usr/bin/env python3
"""Build reconciled, path-free assembly/annotation QC publication tables.

The builder accepts four explicit TSV inputs:

* a metadata/decision table keyed by ``assembly_unit_id`` and containing a
  one-to-one ``qc_sample`` mapping;
* the exact output of :mod:`scripts.qc.basic_stats`, keyed by ``sample``;
* genome-mode and protein-mode exact outputs of
  :mod:`scripts.qc.collect_busco`, also keyed by ``sample``.

No table is published until all identifiers reconcile exactly, every key is
unique, the BUSCO release/dataset/creation date is uniform, and the reported
counts and percentages are arithmetically consistent.  Candidate and excluded
assemblies are deliberately retained.  Private runtime path columns from the
producer tables are never copied to an output.

The complete output directory is assembled in a sibling staging directory and
published with one rename.  An existing nonempty output directory is refused,
so a failed or ambiguous run cannot partly overwrite an accepted result.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import shutil
import sys
import tempfile
from typing import Iterable, Mapping, Sequence


SCRIPT_VERSION = "1.0.0"
PUBLIC_FILE_MODE = 0o644
PUBLIC_DIRECTORY_MODE = 0o755

METADATA_INPUT_COLUMNS = (
    "assembly_unit_id",
    "qc_sample",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "accession",
    "decision_status",
    "publisher_assembly_scope",
    "local_qc_scope",
    "publisher_assembly_provenance",
    "publisher_annotation_provenance",
    "publisher_protein_provenance",
    "decision_reason",
)

# Exact schemas emitted by basic_stats.py and collect_busco.py.  Keeping these
# explicit prevents a similarly named hand-edited table from entering the
# publication join unnoticed.
BASIC_INPUT_COLUMNS = (
    "sample",
    "current_or_alternative",
    "accession",
    "source_url",
    "genome_path",
    "gff_path",
    "protein_path",
    "genome_sequence_count",
    "genome_total_bp",
    "genome_ungapped_bp",
    "genome_n_bp",
    "genome_n_percent",
    "genome_gc_bp",
    "genome_gc_percent",
    "genome_longest_bp",
    "genome_n50_bp",
    "genome_l50",
    "gff_feature_rows",
    "gff_invalid_rows",
    "gff_gene_count",
    "gff_mrna_count",
    "gff_transcript_count",
    "gff_mrna_or_transcript_count",
    "gff_cds_count",
    "gff_exon_count",
    "protein_sequence_count",
    "protein_empty_sequence_count",
    "protein_total_aa",
    "protein_longest_aa",
    "protein_n50_aa",
    "protein_l50",
    "protein_internal_stop_record_count",
    "protein_terminal_stop_record_count",
    "protein_internal_stop_character_count",
    "protein_nonstandard_character_record_count",
    "protein_nonstandard_character_count",
)

BUSCO_INPUT_COLUMNS = (
    "sample",
    "busco_version",
    "dataset",
    "dataset_creation_date",
    "mode",
    "input_path",
    "C_percent",
    "S_percent",
    "D_percent",
    "F_percent",
    "M_percent",
    "n",
    "C_count",
    "S_count",
    "D_count",
    "F_count",
    "M_count",
    "short_summary_path",
)

METADATA_OUTPUT_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "accession",
    "assembly_scope",
    "publisher_assembly_scope",
    "local_qc_scope",
    "decision_status",
    "publisher_assembly_provenance",
    "publisher_annotation_provenance",
    "publisher_protein_provenance",
    "decision_reason",
    "source_table_basename",
    "source_table_sha256",
)

BASIC_PRIVATE_COLUMNS = frozenset(
    {"sample", "genome_path", "gff_path", "protein_path"}
)
BASIC_OUTPUT_COLUMNS = (
    "assembly_unit_id",
    *(column for column in BASIC_INPUT_COLUMNS if column not in BASIC_PRIVATE_COLUMNS),
    "source_table_basename",
    "source_table_sha256",
)

BUSCO_PRIVATE_COLUMNS = frozenset({"sample", "input_path", "short_summary_path"})
BUSCO_OUTPUT_COLUMNS = (
    "assembly_unit_id",
    *(column for column in BUSCO_INPUT_COLUMNS if column not in BUSCO_PRIVATE_COLUMNS),
    "source_table_basename",
    "source_table_sha256",
)

COMBINED_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "accession",
    "decision_status",
    "decision_reason",
    "publisher_assembly_scope",
    "local_qc_scope",
    "publisher_assembly_provenance",
    "publisher_annotation_provenance",
    "publisher_protein_provenance",
    "genome_sequence_count",
    "genome_total_bp",
    "genome_ungapped_bp",
    "genome_n_bp",
    "genome_n_percent",
    "genome_n50_bp",
    "genome_l50",
    "genome_longest_bp",
    "gff_gene_count",
    "gff_mrna_or_transcript_count",
    "protein_sequence_count",
    "genome_busco_version",
    "genome_busco_dataset",
    "genome_busco_dataset_creation_date",
    "genome_busco_mode",
    "genome_busco_C_percent",
    "genome_busco_S_percent",
    "genome_busco_D_percent",
    "genome_busco_F_percent",
    "genome_busco_M_percent",
    "genome_busco_n",
    "genome_busco_C_count",
    "genome_busco_S_count",
    "genome_busco_D_count",
    "genome_busco_F_count",
    "genome_busco_M_count",
    "protein_busco_version",
    "protein_busco_dataset",
    "protein_busco_dataset_creation_date",
    "protein_busco_mode",
    "protein_busco_C_percent",
    "protein_busco_S_percent",
    "protein_busco_D_percent",
    "protein_busco_F_percent",
    "protein_busco_M_percent",
    "protein_busco_n",
    "protein_busco_C_count",
    "protein_busco_S_count",
    "protein_busco_D_count",
    "protein_busco_F_count",
    "protein_busco_M_count",
    "basic_stats_method",
    "basic_stats_source_basename",
    "basic_stats_source_sha256",
    "genome_busco_source_basename",
    "genome_busco_source_sha256",
    "protein_busco_source_basename",
    "protein_busco_source_sha256",
)

OUTPUT_FILENAMES = {
    "metadata": "qc_metadata_public.tsv",
    "basic": "qc_basic_stats_public.tsv",
    "genome_busco": "qc_genome_busco_public.tsv",
    "protein_busco": "qc_protein_busco_public.tsv",
    "combined": "assembly_annotation_qc_supplementary.tsv",
    "validation": "qc_publication_validation.json",
}

DECISION_STATUSES = ("current", "candidate", "excluded")
INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z")
DECIMAL_RE = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
WINDOWS_DRIVE_RE = re.compile(r"^[A-Za-z]:[\\/]")
EMBEDDED_POSIX_PATH_RE = re.compile(
    r"(?:^|[\s=(\"'])/[A-Za-z0-9._-]+(?:/|\Z)"
)


class QCPublicationError(RuntimeError):
    """Raised when QC inputs cannot be published unambiguously."""


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    """Return a lowercase SHA-256 digest without loading the file into memory."""

    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(block_size), b""):
                digest.update(block)
    except OSError as error:
        raise QCPublicationError(f"cannot read input {path}: {error}") from error
    return digest.hexdigest()


def input_provenance(path: Path) -> dict[str, str | int]:
    """Return path-free checksum provenance for one input table."""

    if not path.is_file():
        raise QCPublicationError(f"input is not a regular file: {path}")
    if path.name in {"", ".", ".."} or any(
        character in path.name for character in "\t\r\n"
    ):
        raise QCPublicationError(f"input has an invalid basename: {path}")
    return {
        "basename": path.name,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def read_exact_tsv(
    path: Path, expected_columns: Sequence[str], table_name: str
) -> list[dict[str, str]]:
    """Read a nonempty TSV whose header exactly matches its declared producer."""

    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise QCPublicationError(f"{table_name}: input is not a regular file: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            observed = tuple(reader.fieldnames or ())
            expected = tuple(expected_columns)
            if observed != expected:
                raise QCPublicationError(
                    f"{table_name}: schema mismatch; expected={list(expected)!r}; "
                    f"observed={list(observed)!r}"
                )
            rows: list[dict[str, str]] = []
            for line_number, raw_row in enumerate(reader, start=2):
                if None in raw_row:
                    raise QCPublicationError(
                        f"{table_name}: extra fields at line {line_number}"
                    )
                row = {column: (raw_row[column] or "").strip() for column in expected}
                if not any(row.values()):
                    raise QCPublicationError(
                        f"{table_name}: blank data row at line {line_number}"
                    )
                if any("\n" in value or "\r" in value or "\t" in value for value in row.values()):
                    raise QCPublicationError(
                        f"{table_name}: embedded line break or tab at line {line_number}"
                    )
                rows.append(row)
    except (OSError, csv.Error) as error:
        raise QCPublicationError(f"cannot parse {table_name} {path}: {error}") from error
    if not rows:
        raise QCPublicationError(f"{table_name}: no data rows")
    return rows


def index_unique(
    rows: Iterable[dict[str, str]], key: str, table_name: str
) -> dict[str, dict[str, str]]:
    """Index rows by one required, nonempty, unique key."""

    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row[key]
        if not value:
            raise QCPublicationError(f"{table_name}: empty {key}")
        if value in indexed:
            raise QCPublicationError(f"{table_name}: duplicate {key} {value!r}")
        indexed[value] = row
    return indexed


def parse_integer(value: str, *, field: str, sample: str, positive: bool = False) -> int:
    """Parse a canonical nonnegative integer, optionally requiring positivity."""

    if INTEGER_RE.fullmatch(value) is None:
        raise QCPublicationError(
            f"{sample}: {field} must be a canonical nonnegative integer, found {value!r}"
        )
    parsed = int(value)
    if positive and parsed == 0:
        raise QCPublicationError(f"{sample}: {field} must be greater than zero")
    return parsed


def parse_percentage(value: str, *, field: str, sample: str) -> Decimal:
    """Parse a fixed-point percentage in the closed interval 0--100."""

    if DECIMAL_RE.fullmatch(value) is None:
        raise QCPublicationError(
            f"{sample}: {field} must be a fixed-point percentage, found {value!r}"
        )
    try:
        parsed = Decimal(value)
    except InvalidOperation as error:  # pragma: no cover - regex already guards this
        raise QCPublicationError(f"{sample}: invalid {field} {value!r}") from error
    if not Decimal(0) <= parsed <= Decimal(100):
        raise QCPublicationError(f"{sample}: {field} must be within [0, 100]")
    return parsed


def require_percentage_consistency(
    value: str,
    numerator: int,
    denominator: int,
    *,
    field: str,
    sample: str,
) -> None:
    """Require a reported percentage to match counts at its printed precision."""

    reported = parse_percentage(value, field=field, sample=sample)
    if denominator <= 0:
        raise QCPublicationError(f"{sample}: {field} has a nonpositive denominator")
    expected = Decimal(numerator) * Decimal(100) / Decimal(denominator)
    decimal_places = len(value.partition(".")[2])
    tolerance = Decimal(5).scaleb(-(decimal_places + 1)) + Decimal("1e-12")
    if abs(reported - expected) > tolerance:
        raise QCPublicationError(
            f"{sample}: {field}={value} is inconsistent with "
            f"{numerator}/{denominator} at {decimal_places} decimal places"
        )


def validate_metadata(rows: Sequence[dict[str, str]]) -> None:
    """Validate explicit biological identity, scope, provenance, and decision fields."""

    index_unique(rows, "assembly_unit_id", "metadata")
    index_unique(rows, "qc_sample", "metadata")
    required_nonempty = (
        "biological_species",
        "individual_id",
        "haplotype_or_subgenome",
        "accession",
        "publisher_assembly_scope",
        "local_qc_scope",
        "publisher_assembly_provenance",
        "publisher_annotation_provenance",
        "publisher_protein_provenance",
    )
    for row in rows:
        unit = row["assembly_unit_id"]
        for field in required_nonempty:
            if not row[field]:
                raise QCPublicationError(f"{unit}: metadata field {field} must be nonempty")
        status = row["decision_status"]
        if status not in DECISION_STATUSES:
            raise QCPublicationError(
                f"{unit}: decision_status must be one of {', '.join(DECISION_STATUSES)}"
            )
        if status in {"candidate", "excluded"} and not row["decision_reason"]:
            raise QCPublicationError(
                f"{unit}: {status} rows require a nonempty decision_reason"
            )


def validate_basic_row(row: Mapping[str, str], metadata_row: Mapping[str, str]) -> None:
    """Validate internal arithmetic and annotation-count constraints for one row."""

    sample = row["sample"]
    if row["accession"] != metadata_row["accession"]:
        raise QCPublicationError(
            f"{sample}: basic-statistics accession {row['accession']!r} differs from "
            f"metadata accession {metadata_row['accession']!r}"
        )

    integer_fields = (
        "genome_sequence_count",
        "genome_total_bp",
        "genome_ungapped_bp",
        "genome_n_bp",
        "genome_gc_bp",
        "genome_longest_bp",
        "genome_n50_bp",
        "genome_l50",
        "gff_feature_rows",
        "gff_invalid_rows",
        "gff_gene_count",
        "gff_mrna_count",
        "gff_transcript_count",
        "gff_mrna_or_transcript_count",
        "gff_cds_count",
        "gff_exon_count",
        "protein_sequence_count",
        "protein_empty_sequence_count",
        "protein_total_aa",
        "protein_longest_aa",
        "protein_n50_aa",
        "protein_l50",
        "protein_internal_stop_record_count",
        "protein_terminal_stop_record_count",
        "protein_internal_stop_character_count",
        "protein_nonstandard_character_record_count",
        "protein_nonstandard_character_count",
    )
    values = {
        field: parse_integer(row[field], field=field, sample=sample)
        for field in integer_fields
    }

    sequence_count = values["genome_sequence_count"]
    total_bp = values["genome_total_bp"]
    ungapped_bp = values["genome_ungapped_bp"]
    n_bp = values["genome_n_bp"]
    gc_bp = values["genome_gc_bp"]
    longest_bp = values["genome_longest_bp"]
    n50_bp = values["genome_n50_bp"]
    l50 = values["genome_l50"]
    if sequence_count == 0 or total_bp == 0 or ungapped_bp == 0:
        raise QCPublicationError(f"{sample}: genome FASTA metrics must be nonzero")
    if not (0 < ungapped_bp <= total_bp):
        raise QCPublicationError(f"{sample}: genome_ungapped_bp exceeds genome_total_bp")
    if n_bp + gc_bp > ungapped_bp:
        raise QCPublicationError(f"{sample}: genome N and GC counts exceed ungapped length")
    if not (0 < n50_bp <= longest_bp <= total_bp):
        raise QCPublicationError(f"{sample}: genome N50/longest/total lengths are inconsistent")
    if not 1 <= l50 <= sequence_count:
        raise QCPublicationError(f"{sample}: genome_l50 is outside the sequence-count range")
    require_percentage_consistency(
        row["genome_n_percent"], n_bp, ungapped_bp,
        field="genome_n_percent", sample=sample,
    )
    require_percentage_consistency(
        row["genome_gc_percent"], gc_bp, ungapped_bp,
        field="genome_gc_percent", sample=sample,
    )

    feature_rows = values["gff_feature_rows"]
    for field in (
        "gff_gene_count",
        "gff_mrna_count",
        "gff_transcript_count",
        "gff_cds_count",
        "gff_exon_count",
    ):
        if values[field] > feature_rows:
            raise QCPublicationError(f"{sample}: {field} exceeds gff_feature_rows")
    if sum(
        values[field]
        for field in (
            "gff_gene_count",
            "gff_mrna_count",
            "gff_transcript_count",
            "gff_cds_count",
            "gff_exon_count",
        )
    ) > feature_rows:
        raise QCPublicationError(
            f"{sample}: counted GFF feature types exceed gff_feature_rows"
        )
    if values["gff_mrna_or_transcript_count"] != (
        values["gff_mrna_count"] + values["gff_transcript_count"]
    ):
        raise QCPublicationError(
            f"{sample}: gff_mrna_or_transcript_count is not mRNA plus transcript"
        )

    protein_count = values["protein_sequence_count"]
    empty_count = values["protein_empty_sequence_count"]
    protein_total = values["protein_total_aa"]
    protein_longest = values["protein_longest_aa"]
    protein_n50 = values["protein_n50_aa"]
    protein_l50 = values["protein_l50"]
    if protein_count == 0:
        raise QCPublicationError(f"{sample}: protein_sequence_count must be greater than zero")
    if empty_count > protein_count:
        raise QCPublicationError(f"{sample}: empty protein records exceed protein records")
    if protein_total == 0:
        if (protein_longest, protein_n50, protein_l50, empty_count) != (
            0, 0, 0, protein_count
        ):
            raise QCPublicationError(f"{sample}: empty protein-set metrics are inconsistent")
    else:
        if not (0 < protein_n50 <= protein_longest <= protein_total):
            raise QCPublicationError(
                f"{sample}: protein N50/longest/total lengths are inconsistent"
            )
        if not 1 <= protein_l50 <= protein_count - empty_count:
            raise QCPublicationError(f"{sample}: protein_l50 is outside the usable-record range")
    for field in (
        "protein_internal_stop_record_count",
        "protein_terminal_stop_record_count",
        "protein_nonstandard_character_record_count",
    ):
        if values[field] > protein_count:
            raise QCPublicationError(f"{sample}: {field} exceeds protein_sequence_count")
    for field in (
        "protein_internal_stop_character_count",
        "protein_nonstandard_character_count",
    ):
        if values[field] > protein_total:
            raise QCPublicationError(f"{sample}: {field} exceeds protein_total_aa")


def validate_busco_rows(
    rows: Sequence[dict[str, str]], table_name: str, expected_role: str
) -> tuple[dict[str, dict[str, str]], tuple[str, str, str, str, str]]:
    """Validate BUSCO rows and return normalized rows plus their sole signature."""

    indexed = index_unique(rows, "sample", table_name)
    normalized: dict[str, dict[str, str]] = {}
    signatures: set[tuple[str, str, str, str, str]] = set()
    for sample, row in indexed.items():
        for field in ("busco_version", "dataset", "dataset_creation_date", "mode"):
            if not row[field]:
                raise QCPublicationError(f"{sample}: BUSCO field {field} must be nonempty")
        n = parse_integer(row["n"], field="BUSCO n", sample=sample, positive=True)
        counts = {
            code: parse_integer(
                row[f"{code}_count"], field=f"BUSCO {code}_count", sample=sample
            )
            for code in "SDFM"
        }
        complete = counts["S"] + counts["D"]
        if row["C_count"]:
            reported_complete = parse_integer(
                row["C_count"], field="BUSCO C_count", sample=sample
            )
            if reported_complete != complete:
                raise QCPublicationError(
                    f"{sample}: BUSCO C_count is not S_count plus D_count"
                )
        if complete + counts["F"] + counts["M"] != n:
            raise QCPublicationError(
                f"{sample}: BUSCO S+D+F+M counts do not equal n"
            )
        count_by_code = {"C": complete, **counts}
        for code in "CSDFM":
            require_percentage_consistency(
                row[f"{code}_percent"], count_by_code[code], n,
                field=f"BUSCO {code}_percent", sample=sample,
            )

        normalized_row = dict(row)
        normalized_row["C_count"] = str(complete)
        normalized[sample] = normalized_row
        signatures.add(
            (
                row["busco_version"],
                row["dataset"],
                row["dataset_creation_date"],
                row["mode"],
                row["n"],
            )
        )
    if len(signatures) != 1:
        raise QCPublicationError(
            f"{table_name}: BUSCO version/dataset/date/mode/n is not uniform"
        )
    signature = next(iter(signatures))
    mode = signature[3].lower()
    if expected_role == "protein" and "protein" not in mode:
        raise QCPublicationError(
            f"{table_name}: mode {signature[3]!r} does not identify a protein run"
        )
    if expected_role == "genome" and "protein" in mode:
        raise QCPublicationError(
            f"{table_name}: mode {signature[3]!r} identifies a protein run"
        )
    return normalized, signature


def require_exact_sample_sets(
    expected: set[str], observed: set[str], table_name: str
) -> None:
    """Require a raw producer table to cover every declared QC sample exactly."""

    if observed != expected:
        raise QCPublicationError(
            f"{table_name}: qc_sample set differs from metadata; "
            f"missing={sorted(expected - observed)!r}; extra={sorted(observed - expected)!r}"
        )


def looks_like_local_path(value: str) -> bool:
    """Identify common private local-path forms while allowing public URLs/DOIs."""

    if not value:
        return False
    if value.startswith(("/", "~/", "file://", "\\\\")):
        return True
    if WINDOWS_DRIVE_RE.match(value) or PureWindowsPath(value).is_absolute():
        return True
    return EMBEDDED_POSIX_PATH_RE.search(value) is not None


def reject_path_leaks(groups: Mapping[str, Sequence[Mapping[str, str]]]) -> None:
    """Reject a retained output field that resembles a private runtime path."""

    for table_name, rows in groups.items():
        for row in rows:
            for column, value in row.items():
                if looks_like_local_path(value):
                    raise QCPublicationError(
                        f"{table_name}: retained field {column!r} contains a local path"
                    )


def normalize_tables(
    metadata_rows: Sequence[dict[str, str]],
    basic_rows: Sequence[dict[str, str]],
    genome_rows: Sequence[dict[str, str]],
    protein_rows: Sequence[dict[str, str]],
    provenance: Mapping[str, Mapping[str, str | int]],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, object]]:
    """Validate all joins and return the five exact output tables."""

    validate_metadata(metadata_rows)
    metadata_by_sample = index_unique(metadata_rows, "qc_sample", "metadata")
    expected_samples = set(metadata_by_sample)
    basic_by_sample = index_unique(basic_rows, "sample", "basic_stats")
    genome_by_sample, genome_signature = validate_busco_rows(
        genome_rows, "genome_busco", "genome"
    )
    protein_by_sample, protein_signature = validate_busco_rows(
        protein_rows, "protein_busco", "protein"
    )
    require_exact_sample_sets(expected_samples, set(basic_by_sample), "basic_stats")
    require_exact_sample_sets(expected_samples, set(genome_by_sample), "genome_busco")
    require_exact_sample_sets(expected_samples, set(protein_by_sample), "protein_busco")

    comparable_indices = (0, 1, 2, 4)
    if tuple(genome_signature[index] for index in comparable_indices) != tuple(
        protein_signature[index] for index in comparable_indices
    ):
        raise QCPublicationError(
            "genome and protein BUSCO tables do not share one "
            "version/dataset/creation-date/n signature"
        )
    if genome_signature[3] == protein_signature[3]:
        raise QCPublicationError("genome and protein BUSCO modes must differ")

    normalized_metadata: list[dict[str, str]] = []
    normalized_basic: list[dict[str, str]] = []
    normalized_genome: list[dict[str, str]] = []
    normalized_protein: list[dict[str, str]] = []
    combined: list[dict[str, str]] = []

    for metadata_row in metadata_rows:
        sample = metadata_row["qc_sample"]
        unit = metadata_row["assembly_unit_id"]
        basic = basic_by_sample[sample]
        genome = genome_by_sample[sample]
        protein = protein_by_sample[sample]
        validate_basic_row(basic, metadata_row)

        metadata_output = {
            "assembly_unit_id": unit,
            "biological_species": metadata_row["biological_species"],
            "individual_id": metadata_row["individual_id"],
            "haplotype_or_subgenome": metadata_row["haplotype_or_subgenome"],
            "accession": metadata_row["accession"],
            "assembly_scope": metadata_row["local_qc_scope"],
            "publisher_assembly_scope": metadata_row["publisher_assembly_scope"],
            "local_qc_scope": metadata_row["local_qc_scope"],
            "decision_status": metadata_row["decision_status"],
            "publisher_assembly_provenance": metadata_row[
                "publisher_assembly_provenance"
            ],
            "publisher_annotation_provenance": metadata_row[
                "publisher_annotation_provenance"
            ],
            "publisher_protein_provenance": metadata_row[
                "publisher_protein_provenance"
            ],
            "decision_reason": metadata_row["decision_reason"],
            "source_table_basename": str(provenance["metadata"]["basename"]),
            "source_table_sha256": str(provenance["metadata"]["sha256"]),
        }
        normalized_metadata.append(metadata_output)

        basic_output = {"assembly_unit_id": unit}
        basic_output.update(
            (column, basic[column])
            for column in BASIC_INPUT_COLUMNS
            if column not in BASIC_PRIVATE_COLUMNS
        )
        basic_output.update(
            {
                "source_table_basename": str(provenance["basic"]["basename"]),
                "source_table_sha256": str(provenance["basic"]["sha256"]),
            }
        )
        normalized_basic.append(basic_output)

        genome_output = {"assembly_unit_id": unit}
        genome_output.update(
            (column, genome[column])
            for column in BUSCO_INPUT_COLUMNS
            if column not in BUSCO_PRIVATE_COLUMNS
        )
        genome_output.update(
            {
                "source_table_basename": str(provenance["genome_busco"]["basename"]),
                "source_table_sha256": str(provenance["genome_busco"]["sha256"]),
            }
        )
        normalized_genome.append(genome_output)

        protein_output = {"assembly_unit_id": unit}
        protein_output.update(
            (column, protein[column])
            for column in BUSCO_INPUT_COLUMNS
            if column not in BUSCO_PRIVATE_COLUMNS
        )
        protein_output.update(
            {
                "source_table_basename": str(provenance["protein_busco"]["basename"]),
                "source_table_sha256": str(provenance["protein_busco"]["sha256"]),
            }
        )
        normalized_protein.append(protein_output)

        combined_row = {
            "assembly_unit_id": unit,
            "biological_species": metadata_row["biological_species"],
            "individual_id": metadata_row["individual_id"],
            "haplotype_or_subgenome": metadata_row["haplotype_or_subgenome"],
            "accession": metadata_row["accession"],
            "decision_status": metadata_row["decision_status"],
            "decision_reason": metadata_row["decision_reason"],
            "publisher_assembly_scope": metadata_row["publisher_assembly_scope"],
            "local_qc_scope": metadata_row["local_qc_scope"],
            "publisher_assembly_provenance": metadata_row[
                "publisher_assembly_provenance"
            ],
            "publisher_annotation_provenance": metadata_row[
                "publisher_annotation_provenance"
            ],
            "publisher_protein_provenance": metadata_row[
                "publisher_protein_provenance"
            ],
            "genome_sequence_count": basic["genome_sequence_count"],
            "genome_total_bp": basic["genome_total_bp"],
            "genome_ungapped_bp": basic["genome_ungapped_bp"],
            "genome_n_bp": basic["genome_n_bp"],
            "genome_n_percent": basic["genome_n_percent"],
            "genome_n50_bp": basic["genome_n50_bp"],
            "genome_l50": basic["genome_l50"],
            "genome_longest_bp": basic["genome_longest_bp"],
            "gff_gene_count": basic["gff_gene_count"],
            "gff_mrna_or_transcript_count": basic["gff_mrna_or_transcript_count"],
            "protein_sequence_count": basic["protein_sequence_count"],
            "genome_busco_version": genome["busco_version"],
            "genome_busco_dataset": genome["dataset"],
            "genome_busco_dataset_creation_date": genome["dataset_creation_date"],
            "genome_busco_mode": genome["mode"],
            "protein_busco_version": protein["busco_version"],
            "protein_busco_dataset": protein["dataset"],
            "protein_busco_dataset_creation_date": protein["dataset_creation_date"],
            "protein_busco_mode": protein["mode"],
            "basic_stats_method": "basic_stats.py",
            "basic_stats_source_basename": str(provenance["basic"]["basename"]),
            "basic_stats_source_sha256": str(provenance["basic"]["sha256"]),
            "genome_busco_source_basename": str(
                provenance["genome_busco"]["basename"]
            ),
            "genome_busco_source_sha256": str(provenance["genome_busco"]["sha256"]),
            "protein_busco_source_basename": str(
                provenance["protein_busco"]["basename"]
            ),
            "protein_busco_source_sha256": str(
                provenance["protein_busco"]["sha256"]
            ),
        }
        for prefix, busco in (("genome", genome), ("protein", protein)):
            for code in "CSDFM":
                combined_row[f"{prefix}_busco_{code}_percent"] = busco[
                    f"{code}_percent"
                ]
            combined_row[f"{prefix}_busco_n"] = busco["n"]
            for code in "CSDFM":
                combined_row[f"{prefix}_busco_{code}_count"] = busco[f"{code}_count"]
        combined.append(combined_row)

    groups = {
        "metadata": normalized_metadata,
        "basic": normalized_basic,
        "genome_busco": normalized_genome,
        "protein_busco": normalized_protein,
        "combined": combined,
    }
    reject_path_leaks(groups)
    status_counts = {
        status: sum(row["decision_status"] == status for row in normalized_metadata)
        for status in DECISION_STATUSES
    }
    validation: dict[str, object] = {
        "schema_version": "1.0",
        "builder": "build_qc_publication_tables.py",
        "builder_version": SCRIPT_VERSION,
        "status": "pass",
        "assembly_unit_count": len(metadata_rows),
        "decision_status_counts": status_counts,
        "busco_signature": {
            "version": genome_signature[0],
            "dataset": genome_signature[1],
            "dataset_creation_date": genome_signature[2],
            "n": genome_signature[4],
            "genome_mode": genome_signature[3],
            "protein_mode": protein_signature[3],
        },
        "checks": {
            "assembly_unit_ids_unique": "pass",
            "qc_samples_unique": "pass",
            "producer_sample_sets_exactly_match_metadata": "pass",
            "all_rows_retained_in_metadata_order": "pass",
            "candidate_and_excluded_rows_retained": "pass",
            "basic_statistics_arithmetic": "pass",
            "busco_version_dataset_date_and_n_uniform": "pass",
            "busco_counts_and_percentages_consistent": "pass",
            "private_runtime_columns_removed": "pass",
            "retained_fields_path_free": "pass",
        },
        "inputs": {
            role: dict(entry) for role, entry in provenance.items()
        },
    }
    return groups, validation


def tsv_bytes(columns: Sequence[str], rows: Sequence[Mapping[str, str]]) -> bytes:
    """Serialize exact-schema rows deterministically as UTF-8 TSV bytes."""

    import io

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(columns),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def write_staged_file(path: Path, content: bytes) -> None:
    """Write and fsync one file inside a private staging directory."""

    with path.open("wb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(path, PUBLIC_FILE_MODE)


def prepare_destination(output: Path) -> bool:
    """Validate an output location and report whether an empty directory exists."""

    if output.is_symlink():
        raise QCPublicationError(f"refusing symlink output directory: {output}")
    if not output.exists():
        return False
    if not output.is_dir():
        raise QCPublicationError(f"output exists and is not a directory: {output}")
    if any(output.iterdir()):
        raise QCPublicationError(f"refusing nonempty output directory: {output}")
    return True


def publish_bundle(
    output_dir: Path,
    groups: Mapping[str, Sequence[Mapping[str, str]]],
    validation: Mapping[str, object],
) -> Path:
    """Atomically publish all normalized TSVs and their validation report."""

    output = Path(output_dir).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    existing_empty = prepare_destination(output)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        os.chmod(stage, PUBLIC_DIRECTORY_MODE)
        schemas = {
            "metadata": METADATA_OUTPUT_COLUMNS,
            "basic": BASIC_OUTPUT_COLUMNS,
            "genome_busco": BUSCO_OUTPUT_COLUMNS,
            "protein_busco": BUSCO_OUTPUT_COLUMNS,
            "combined": COMBINED_COLUMNS,
        }
        for role, columns in schemas.items():
            write_staged_file(
                stage / OUTPUT_FILENAMES[role],
                tsv_bytes(columns, groups[role]),
            )
        validation_content = (
            json.dumps(validation, ensure_ascii=False, allow_nan=False, indent=2, sort_keys=True)
            + "\n"
        ).encode("utf-8")
        write_staged_file(stage / OUTPUT_FILENAMES["validation"], validation_content)
        directory_fd = os.open(stage, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        if existing_empty:
            output.rmdir()
        os.replace(stage, output)
        parent_fd = os.open(output.parent, os.O_RDONLY)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--basic-stats", required=True, type=Path)
    parser.add_argument("--genome-busco", required=True, type=Path)
    parser.add_argument("--protein-busco", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    """Validate all four inputs and publish a complete path-free output bundle."""

    paths = {
        "metadata": Path(args.metadata).expanduser().resolve(),
        "basic": Path(args.basic_stats).expanduser().resolve(),
        "genome_busco": Path(args.genome_busco).expanduser().resolve(),
        "protein_busco": Path(args.protein_busco).expanduser().resolve(),
    }
    provenance = {role: input_provenance(path) for role, path in paths.items()}
    metadata_rows = read_exact_tsv(
        paths["metadata"], METADATA_INPUT_COLUMNS, "metadata"
    )
    basic_rows = read_exact_tsv(paths["basic"], BASIC_INPUT_COLUMNS, "basic_stats")
    genome_rows = read_exact_tsv(
        paths["genome_busco"], BUSCO_INPUT_COLUMNS, "genome_busco"
    )
    protein_rows = read_exact_tsv(
        paths["protein_busco"], BUSCO_INPUT_COLUMNS, "protein_busco"
    )
    groups, validation = normalize_tables(
        metadata_rows, basic_rows, genome_rows, protein_rows, provenance
    )
    output = publish_bundle(args.output_dir, groups, validation)
    return {
        "status": "complete",
        "assembly_unit_count": len(groups["metadata"]),
        "output_directory_name": output.name,
        "outputs": OUTPUT_FILENAMES,
    }


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = run(args)
    except (QCPublicationError, OSError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
