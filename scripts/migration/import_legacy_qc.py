#!/usr/bin/env python3
"""Import completed legacy assembly-QC summaries as public, path-free TSVs.

The importer accepts one or more files for each of three table types: basic
assembly/annotation statistics, genome BUSCO summaries, and protein BUSCO
summaries.  It validates every input before writing any output.  In
particular, legacy sample labels must be mapped explicitly, rows may not be
duplicated across split input files, input schemas must match the known
producer schemas exactly, and all BUSCO tables must use one compatible BUSCO
version and lineage-dataset signature.

Absolute runtime path columns are omitted from the public tables.  Each output
row instead records the basename and SHA256 checksum of the exact source TSV
from which it was imported.  This preserves traceability without publishing
private filesystem locations or copying large biological data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import sys
import tempfile
from typing import Iterable, Mapping, Sequence


MAPPING_COLUMNS = (
    "legacy_sample",
    "assembly_unit_id",
)

BASIC_STATS_COLUMNS = (
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

BUSCO_COLUMNS = (
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

BASIC_PRIVATE_COLUMNS = frozenset({"sample", "genome_path", "gff_path", "protein_path"})
BUSCO_PRIVATE_COLUMNS = frozenset({"sample", "input_path", "short_summary_path"})
PROVENANCE_COLUMNS = ("source_basename", "source_sha256")

PUBLIC_BASIC_COLUMNS = (
    "assembly_unit_id",
    "legacy_sample",
    *(column for column in BASIC_STATS_COLUMNS if column not in BASIC_PRIVATE_COLUMNS),
    *PROVENANCE_COLUMNS,
)
PUBLIC_BUSCO_COLUMNS = (
    "assembly_unit_id",
    "legacy_sample",
    *(column for column in BUSCO_COLUMNS if column not in BUSCO_PRIVATE_COLUMNS),
    *PROVENANCE_COLUMNS,
)

OUTPUT_FILENAMES = {
    "basic_stats": "legacy_qc_basic_stats_public.tsv",
    "genome_busco": "legacy_qc_genome_busco_public.tsv",
    "protein_busco": "legacy_qc_protein_busco_public.tsv",
}


class LegacyQCImportError(RuntimeError):
    """Raised when legacy QC tables cannot be imported unambiguously."""


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    """Return the lowercase SHA256 digest of a file."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(block_size), b""):
                digest.update(block)
    except OSError as error:
        raise LegacyQCImportError(f"cannot read {path}: {error}") from error
    return digest.hexdigest()


def read_exact_tsv(path: Path, expected_columns: Sequence[str], table_name: str) -> list[dict[str, str]]:
    """Read one TSV after requiring the exact producer schema and nonempty rows."""
    path = Path(path)
    if not path.is_file():
        raise LegacyQCImportError(f"{table_name}: input is not a regular file: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            observed = tuple(reader.fieldnames or ())
            expected = tuple(expected_columns)
            if observed != expected:
                raise LegacyQCImportError(
                    f"{table_name}: schema mismatch in {path}; "
                    f"expected={list(expected)!r}; observed={list(observed)!r}"
                )
            rows: list[dict[str, str]] = []
            for line_number, raw_row in enumerate(reader, start=2):
                if None in raw_row:
                    raise LegacyQCImportError(
                        f"{table_name}: extra fields at {path}:{line_number}"
                    )
                row = {column: (raw_row[column] or "").strip() for column in expected}
                if not any(row.values()):
                    raise LegacyQCImportError(
                        f"{table_name}: blank data row at {path}:{line_number}"
                    )
                rows.append(row)
    except (OSError, csv.Error) as error:
        raise LegacyQCImportError(f"cannot parse {table_name} input {path}: {error}") from error
    if not rows:
        raise LegacyQCImportError(f"{table_name}: no data rows in {path}")
    return rows


def load_mapping(path: Path) -> dict[str, str]:
    """Load the public one-to-one legacy-label mapping."""
    rows = read_exact_tsv(path, MAPPING_COLUMNS, "mapping")
    mapping: dict[str, str] = {}
    used_units: dict[str, str] = {}
    for row in rows:
        legacy_sample = row["legacy_sample"]
        assembly_unit_id = row["assembly_unit_id"]
        if not legacy_sample or not assembly_unit_id:
            raise LegacyQCImportError("mapping: legacy_sample and assembly_unit_id must be nonempty")
        if legacy_sample in mapping:
            raise LegacyQCImportError(f"mapping: duplicate legacy_sample {legacy_sample!r}")
        if assembly_unit_id in used_units:
            raise LegacyQCImportError(
                "mapping: duplicate assembly_unit_id "
                f"{assembly_unit_id!r} for {used_units[assembly_unit_id]!r} and {legacy_sample!r}"
            )
        mapping[legacy_sample] = assembly_unit_id
        used_units[assembly_unit_id] = legacy_sample
    return mapping


def path_provenance(path: Path) -> dict[str, str]:
    """Return public provenance for an input without exposing its directory."""
    return {
        "source_basename": path.name,
        "source_sha256": sha256_file(path),
    }


def import_group(
    paths: Iterable[Path],
    *,
    expected_columns: Sequence[str],
    private_columns: frozenset[str],
    table_name: str,
    mapping: Mapping[str, str],
) -> list[dict[str, str]]:
    """Read, map, and sanitize all files belonging to one QC table type."""
    imported: list[dict[str, str]] = []
    seen_samples: dict[str, str] = {}
    for original_path in paths:
        path = Path(original_path)
        provenance = path_provenance(path)
        rows = read_exact_tsv(path, expected_columns, table_name)
        for row in rows:
            legacy_sample = row["sample"]
            if not legacy_sample:
                raise LegacyQCImportError(f"{table_name}: empty sample in {path}")
            if legacy_sample not in mapping:
                raise LegacyQCImportError(
                    f"{table_name}: unmapped legacy sample {legacy_sample!r} in {path}"
                )
            if legacy_sample in seen_samples:
                raise LegacyQCImportError(
                    f"{table_name}: duplicate sample {legacy_sample!r} in "
                    f"{seen_samples[legacy_sample]} and {path}"
                )
            seen_samples[legacy_sample] = str(path)
            public_row = {
                "assembly_unit_id": mapping[legacy_sample],
                "legacy_sample": legacy_sample,
            }
            public_row.update(
                (column, row[column])
                for column in expected_columns
                if column not in private_columns
            )
            public_row.update(provenance)
            imported.append(public_row)
    return sorted(imported, key=lambda row: (row["assembly_unit_id"], row["legacy_sample"]))


def require_same_sample_sets(groups: Mapping[str, Sequence[Mapping[str, str]]]) -> None:
    """Require the three imported table types to cover exactly the same units."""
    sample_sets = {
        name: {row["legacy_sample"] for row in rows}
        for name, rows in groups.items()
    }
    names = tuple(sample_sets)
    reference_name = names[0]
    reference = sample_sets[reference_name]
    for name in names[1:]:
        observed = sample_sets[name]
        if observed != reference:
            raise LegacyQCImportError(
                "QC sample sets differ: "
                f"{reference_name}_only={sorted(reference - observed)!r}; "
                f"{name}_only={sorted(observed - reference)!r}"
            )


def busco_signature(rows: Sequence[Mapping[str, str]], table_name: str) -> tuple[str, ...]:
    """Return the sole BUSCO run signature, rejecting mixed versions or datasets."""
    signature_columns = ("busco_version", "dataset", "dataset_creation_date", "mode", "n")
    signatures = {
        tuple(row[column] for column in signature_columns)
        for row in rows
    }
    if any(not value for signature in signatures for value in signature):
        raise LegacyQCImportError(
            f"{table_name}: BUSCO version/dataset/mode/n signature contains an empty value"
        )
    if len(signatures) != 1:
        raise LegacyQCImportError(
            f"{table_name}: BUSCO version or dataset signatures differ: {sorted(signatures)!r}"
        )
    return next(iter(signatures))


def validate_busco_compatibility(
    genome_rows: Sequence[Mapping[str, str]],
    protein_rows: Sequence[Mapping[str, str]],
) -> None:
    """Require genome and protein BUSCO results to share version and lineage data."""
    genome = busco_signature(genome_rows, "genome_busco")
    protein = busco_signature(protein_rows, "protein_busco")
    # Modes intentionally differ, so compare version, dataset, creation date,
    # and lineage BUSCO count while leaving the mode at index 3 out.
    comparable_indices = (0, 1, 2, 4)
    if tuple(genome[index] for index in comparable_indices) != tuple(
        protein[index] for index in comparable_indices
    ):
        raise LegacyQCImportError(
            "genome_busco and protein_busco use incompatible BUSCO "
            f"version/dataset signatures: genome={genome!r}; protein={protein!r}"
        )
    if genome[3] == protein[3]:
        raise LegacyQCImportError(
            "genome_busco and protein_busco unexpectedly use the same BUSCO mode "
            f"{genome[3]!r}; check that the inputs were assigned to the correct options"
        )


WINDOWS_DRIVE_RE = re.compile(r"^[A-Za-z]:[\\/]")


def looks_like_absolute_path(value: str) -> bool:
    """Return true for common POSIX, home-relative, file-URL, or Windows paths."""
    if not value:
        return False
    return (
        value.startswith(("/", "~/", "file://", "\\\\"))
        or bool(WINDOWS_DRIVE_RE.match(value))
        or PureWindowsPath(value).is_absolute()
    )


def reject_public_path_leaks(groups: Mapping[str, Sequence[Mapping[str, str]]]) -> None:
    """Fail if a retained public field still looks like a local absolute path."""
    for table_name, rows in groups.items():
        for row in rows:
            for column, value in row.items():
                if looks_like_absolute_path(value):
                    raise LegacyQCImportError(
                        f"{table_name}: retained field {column!r} for "
                        f"{row['legacy_sample']!r} contains an absolute path"
                    )


def write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    """Write a complete TSV to a temporary file and atomically replace its target."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(columns),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.chmod(temporary, 0o644)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path, help="Public legacy-label map TSV")
    parser.add_argument(
        "--basic-stats",
        required=True,
        action="append",
        type=Path,
        help="Basic-statistics TSV; repeat for split cohorts",
    )
    parser.add_argument(
        "--genome-busco",
        required=True,
        action="append",
        type=Path,
        help="Genome-BUSCO TSV; repeat for split cohorts",
    )
    parser.add_argument(
        "--protein-busco",
        required=True,
        action="append",
        type=Path,
        help="Protein-BUSCO TSV; repeat for split cohorts",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for public TSVs")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    """Validate all inputs, then emit the three public summaries."""
    mapping = load_mapping(args.mapping)
    groups = {
        "basic_stats": import_group(
            args.basic_stats,
            expected_columns=BASIC_STATS_COLUMNS,
            private_columns=BASIC_PRIVATE_COLUMNS,
            table_name="basic_stats",
            mapping=mapping,
        ),
        "genome_busco": import_group(
            args.genome_busco,
            expected_columns=BUSCO_COLUMNS,
            private_columns=BUSCO_PRIVATE_COLUMNS,
            table_name="genome_busco",
            mapping=mapping,
        ),
        "protein_busco": import_group(
            args.protein_busco,
            expected_columns=BUSCO_COLUMNS,
            private_columns=BUSCO_PRIVATE_COLUMNS,
            table_name="protein_busco",
            mapping=mapping,
        ),
    }
    require_same_sample_sets(groups)
    validate_busco_compatibility(groups["genome_busco"], groups["protein_busco"])
    reject_public_path_leaks(groups)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / OUTPUT_FILENAMES["basic_stats"],
        PUBLIC_BASIC_COLUMNS,
        groups["basic_stats"],
    )
    write_tsv(
        output_dir / OUTPUT_FILENAMES["genome_busco"],
        PUBLIC_BUSCO_COLUMNS,
        groups["genome_busco"],
    )
    write_tsv(
        output_dir / OUTPUT_FILENAMES["protein_busco"],
        PUBLIC_BUSCO_COLUMNS,
        groups["protein_busco"],
    )
    return {
        "status": "complete",
        "sample_count": len(groups["basic_stats"]),
        "output_directory_name": output_dir.resolve().name,
        "outputs": OUTPUT_FILENAMES,
    }


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = run(args)
    except (LegacyQCImportError, OSError, csv.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
