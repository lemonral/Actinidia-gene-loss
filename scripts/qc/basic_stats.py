#!/usr/bin/env python3
"""Summarize assembly, annotation, and protein-set QC metrics.

This utility intentionally uses only the Python standard library.  It streams
plain-text or gzip-compressed FASTA and GFF3 files, so it can be used directly
on large genome assets without loading complete sequences into memory.

Input manifest
--------------
The input is a tab-separated file with these columns:

    sample  current_or_alternative  accession  genome  gff  protein  source_url

``sample`` and the three asset paths are required.  Metadata fields may be
blank.  Relative asset paths are resolved relative to the manifest file.

Example
-------
python scripts/assembly_qc/basic_stats.py \
  --manifest config/assembly_qc_manifest.tsv \
  --output results/assembly_qc/basic_stats.tsv

Metric conventions
------------------
* FASTA sequence characters are streamed after whitespace removal.
* Genome ``total_bp``, ``longest_bp``, N50, and L50 use full record lengths;
  ``ungapped_bp`` removes ``-`` and ``.`` characters.
* ``N%`` and ``GC%`` use ungapped genome length as the denominator.  GC%
  therefore includes ambiguous non-N characters in the denominator.
* Protein totals and protein length statistics use ungapped sequence lengths.
  A protein record consisting only of gaps is treated as empty.
* Protein stop-codon metrics ignore ``-`` and ``.`` gaps.  One or more ``*``
  characters at the end of the ungapped record are terminal stops; any other
  ``*`` characters are internal stops.  Terminal stops remain included in the
  existing protein length statistics.
* The accepted protein alphabet is the 20 standard amino-acid symbols plus
  ``B``, ``J``, ``O``, ``U``, ``X``, and ``Z``.  Matching is case-insensitive;
  ``*`` is handled as a stop, and ``-``/``.`` gaps are ignored.  Every other
  non-whitespace character is counted as nonstandard.
* GFF3 feature counts are case-insensitive exact matches of feature column 3.
  ``mRNA_or_transcript`` is the sum of separate ``mRNA`` and ``transcript``
  counts; it does not attempt to de-duplicate feature IDs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


REQUIRED_MANIFEST_COLUMNS = (
    "sample",
    "current_or_alternative",
    "accession",
    "genome",
    "gff",
    "protein",
    "source_url",
)

OUTPUT_COLUMNS = (
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

PROTEIN_ALLOWED_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBJOUXZ"
PROTEIN_ACCEPTED_DELETE_TABLE = str.maketrans(
    "",
    "",
    PROTEIN_ALLOWED_ALPHABET + "*-.",
)


class QCInputError(RuntimeError):
    """Raised for a malformed manifest or missing/unreadable input."""


@dataclass(frozen=True)
class FastaSummary:
    sequence_count: int
    total_length: int
    ungapped_length: int
    n_count: int
    gc_count: int
    longest_raw_length: int
    longest_ungapped_length: int
    raw_n50: int
    raw_l50: int
    ungapped_n50: int
    ungapped_l50: int
    empty_ungapped_records: int
    internal_stop_records: int
    terminal_stop_records: int
    internal_stop_characters: int
    nonstandard_character_records: int
    nonstandard_characters: int


@dataclass(frozen=True)
class GffSummary:
    feature_rows: int
    invalid_rows: int
    gene_count: int
    mrna_count: int
    transcript_count: int
    cds_count: int
    exon_count: int


def open_text(path: Path) -> TextIO:
    """Open plain-text or gzip-compressed text, preserving line streaming."""
    try:
        if path.suffix.lower() == ".gz":
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return path.open("rt", encoding="utf-8", errors="replace")
    except OSError as error:
        raise QCInputError(f"Cannot open {path}: {error}") from error


def n50_l50(lengths: list[int]) -> tuple[int, int]:
    """Return N50/L50 from non-zero sequence lengths, or (0, 0) if empty."""
    usable = sorted((length for length in lengths if length > 0), reverse=True)
    total = sum(usable)
    if total == 0:
        return 0, 0

    cumulative = 0
    threshold = total / 2
    for index, length in enumerate(usable, start=1):
        cumulative += length
        if cumulative >= threshold:
            return length, index
    raise AssertionError("N50 calculation did not reach the total sequence length")


def finish_record(
    raw_lengths: list[int],
    ungapped_lengths: list[int],
    raw_length: int,
    ungapped_length: int,
) -> None:
    """Store the two length conventions for one completed FASTA record."""
    raw_lengths.append(raw_length)
    ungapped_lengths.append(ungapped_length)


def summarize_fasta(path: Path, *, protein_integrity: bool = False) -> FastaSummary:
    """Stream a FASTA file and calculate composition, length, and integrity metrics."""
    raw_lengths: list[int] = []
    ungapped_lengths: list[int] = []
    current_record = False
    raw_length = 0
    ungapped_length = 0
    n_count = 0
    gc_count = 0
    record_stop_count = 0
    record_trailing_stop_count = 0
    record_nonstandard_count = 0
    internal_stop_records = 0
    terminal_stop_records = 0
    internal_stop_characters = 0
    nonstandard_character_records = 0
    nonstandard_characters = 0

    def finish_current_record() -> None:
        """Commit the current record without retaining its sequence in memory."""
        nonlocal internal_stop_records
        nonlocal terminal_stop_records
        nonlocal internal_stop_characters
        nonlocal nonstandard_character_records
        nonlocal nonstandard_characters

        finish_record(
            raw_lengths,
            ungapped_lengths,
            raw_length,
            ungapped_length,
        )
        if not protein_integrity:
            return
        record_internal_stops = record_stop_count - record_trailing_stop_count
        internal_stop_records += int(record_internal_stops > 0)
        terminal_stop_records += int(record_trailing_stop_count > 0)
        internal_stop_characters += record_internal_stops
        nonstandard_character_records += int(record_nonstandard_count > 0)
        nonstandard_characters += record_nonstandard_count

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current_record:
                    finish_current_record()
                current_record = True
                raw_length = 0
                ungapped_length = 0
                record_stop_count = 0
                record_trailing_stop_count = 0
                record_nonstandard_count = 0
                continue

            if not current_record:
                raise QCInputError(f"{path}: sequence data before the first FASTA header at line {line_number}")

            sequence = "".join(stripped.split())
            raw_length += len(sequence)
            upper = sequence.upper()
            # Built-in string counts run in C and keep this streaming pass fast
            # enough for multi-gigabase cohorts.  Avoid a Python-level loop over
            # every nucleotide/residue.
            ungapped_length += len(upper) - upper.count("-") - upper.count(".")
            n_count += upper.count("N")
            gc_count += upper.count("G") + upper.count("C")
            if protein_integrity:
                record_stop_count += upper.count("*")
                meaningful = upper.replace("-", "").replace(".", "")
                if meaningful:
                    trailing_stops = len(meaningful) - len(meaningful.rstrip("*"))
                    if trailing_stops == len(meaningful):
                        record_trailing_stop_count += trailing_stops
                    else:
                        record_trailing_stop_count = trailing_stops
                record_nonstandard_count += len(
                    upper.translate(PROTEIN_ACCEPTED_DELETE_TABLE)
                )

    if not current_record:
        raise QCInputError(f"{path}: no FASTA records found")

    finish_current_record()
    raw_n50, raw_l50 = n50_l50(raw_lengths)
    ungapped_n50, ungapped_l50 = n50_l50(ungapped_lengths)
    return FastaSummary(
        sequence_count=len(raw_lengths),
        total_length=sum(raw_lengths),
        ungapped_length=sum(ungapped_lengths),
        n_count=n_count,
        gc_count=gc_count,
        longest_raw_length=max(raw_lengths, default=0),
        longest_ungapped_length=max(ungapped_lengths, default=0),
        raw_n50=raw_n50,
        raw_l50=raw_l50,
        ungapped_n50=ungapped_n50,
        ungapped_l50=ungapped_l50,
        empty_ungapped_records=sum(length == 0 for length in ungapped_lengths),
        internal_stop_records=internal_stop_records,
        terminal_stop_records=terminal_stop_records,
        internal_stop_characters=internal_stop_characters,
        nonstandard_character_records=nonstandard_character_records,
        nonstandard_characters=nonstandard_characters,
    )


def summarize_gff(path: Path) -> GffSummary:
    """Stream GFF3 rows and count core annotation feature types."""
    feature_rows = 0
    invalid_rows = 0
    counts = {"gene": 0, "mrna": 0, "transcript": 0, "cds": 0, "exon": 0}

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                invalid_rows += 1
                continue
            feature_rows += 1
            feature = fields[2].strip().lower()
            if feature in counts:
                counts[feature] += 1

    return GffSummary(
        feature_rows=feature_rows,
        invalid_rows=invalid_rows,
        gene_count=counts["gene"],
        mrna_count=counts["mrna"],
        transcript_count=counts["transcript"],
        cds_count=counts["cds"],
        exon_count=counts["exon"],
    )


def meaningful_manifest_lines(handle: TextIO) -> Iterator[str]:
    """Ignore blank and comment lines while retaining TSV header/data rows."""
    for line in handle:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        yield line


def resolve_asset(raw_path: str, manifest_path: Path, field: str, sample: str) -> Path:
    """Resolve and validate a required asset path from one manifest row."""
    if not raw_path:
        raise QCInputError(f"{manifest_path}: sample {sample!r} has an empty {field} field")
    candidate = Path(raw_path).expanduser()
    if not candidate.is_absolute():
        candidate = manifest_path.parent / candidate
    resolved = candidate.resolve()
    if not resolved.is_file():
        raise QCInputError(f"{manifest_path}: sample {sample!r} {field} is not a readable file: {resolved}")
    return resolved


def read_manifest(path: Path) -> list[dict[str, str | Path]]:
    """Read and validate a manifest, resolving asset paths relative to it."""
    if not path.is_file():
        raise QCInputError(f"Manifest does not exist: {path}")

    with open_text(path) as handle:
        reader = csv.DictReader(meaningful_manifest_lines(handle), delimiter="\t")
        if reader.fieldnames is None:
            raise QCInputError(f"{path}: missing a TSV header")
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        missing = [field for field in REQUIRED_MANIFEST_COLUMNS if field not in reader.fieldnames]
        if missing:
            raise QCInputError(f"{path}: missing required columns: {', '.join(missing)}")

        rows: list[dict[str, str | Path]] = []
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise QCInputError(f"{path}: line {line_number} has more fields than the header")
            row = {key: (value or "").strip() for key, value in raw_row.items()}
            if not any(row.values()):
                continue
            sample = row["sample"]
            if not sample:
                raise QCInputError(f"{path}: line {line_number} has an empty sample")
            row["genome"] = resolve_asset(row["genome"], path, "genome", sample)
            row["gff"] = resolve_asset(row["gff"], path, "gff", sample)
            row["protein"] = resolve_asset(row["protein"], path, "protein", sample)
            rows.append(row)

    if not rows:
        raise QCInputError(f"{path}: no data rows found")
    return rows


def percentage(numerator: int, denominator: int) -> str:
    """Format a percentage without hiding an undefined zero-length denominator."""
    if denominator == 0:
        return "NA"
    return f"{100 * numerator / denominator:.6f}"


def build_output_row(manifest_row: dict[str, str | Path]) -> dict[str, str | int]:
    """Calculate all requested metrics for one sample/assembly row."""
    genome_path = manifest_row["genome"]
    gff_path = manifest_row["gff"]
    protein_path = manifest_row["protein"]
    assert isinstance(genome_path, Path)
    assert isinstance(gff_path, Path)
    assert isinstance(protein_path, Path)

    genome = summarize_fasta(genome_path)
    gff = summarize_gff(gff_path)
    protein = summarize_fasta(protein_path, protein_integrity=True)
    return {
        "sample": str(manifest_row["sample"]),
        "current_or_alternative": str(manifest_row["current_or_alternative"]),
        "accession": str(manifest_row["accession"]),
        "source_url": str(manifest_row["source_url"]),
        "genome_path": str(genome_path),
        "gff_path": str(gff_path),
        "protein_path": str(protein_path),
        "genome_sequence_count": genome.sequence_count,
        "genome_total_bp": genome.total_length,
        "genome_ungapped_bp": genome.ungapped_length,
        "genome_n_bp": genome.n_count,
        "genome_n_percent": percentage(genome.n_count, genome.ungapped_length),
        "genome_gc_bp": genome.gc_count,
        "genome_gc_percent": percentage(genome.gc_count, genome.ungapped_length),
        "genome_longest_bp": genome.longest_raw_length,
        "genome_n50_bp": genome.raw_n50,
        "genome_l50": genome.raw_l50,
        "gff_feature_rows": gff.feature_rows,
        "gff_invalid_rows": gff.invalid_rows,
        "gff_gene_count": gff.gene_count,
        "gff_mrna_count": gff.mrna_count,
        "gff_transcript_count": gff.transcript_count,
        "gff_mrna_or_transcript_count": gff.mrna_count + gff.transcript_count,
        "gff_cds_count": gff.cds_count,
        "gff_exon_count": gff.exon_count,
        "protein_sequence_count": protein.sequence_count,
        "protein_empty_sequence_count": protein.empty_ungapped_records,
        "protein_total_aa": protein.ungapped_length,
        "protein_longest_aa": protein.longest_ungapped_length,
        "protein_n50_aa": protein.ungapped_n50,
        "protein_l50": protein.ungapped_l50,
        "protein_internal_stop_record_count": protein.internal_stop_records,
        "protein_terminal_stop_record_count": protein.terminal_stop_records,
        "protein_internal_stop_character_count": protein.internal_stop_characters,
        "protein_nonstandard_character_record_count": protein.nonstandard_character_records,
        "protein_nonstandard_character_count": protein.nonstandard_characters,
    }


def write_rows(rows: list[dict[str, str | int]], output: str) -> None:
    """Write a deterministic TSV to stdout or an explicit output path."""
    if output == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return

    output_path = Path(output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Stream basic assembly, annotation, and protein-set QC metrics from a TSV manifest."
    )
    parser.add_argument("--manifest", required=True, help="TSV with sample/current_or_alternative/accession/genome/gff/protein/source_url")
    parser.add_argument("--output", default="-", help="Output TSV path, or '-' for stdout (default: '-')")
    parser.add_argument("--version", action="version", version="basic_stats 1.1.0")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest_path = Path(args.manifest).expanduser().resolve()
        manifest_rows = read_manifest(manifest_path)
        output_rows = [build_output_row(row) for row in manifest_rows]
        write_rows(output_rows, args.output)
    except QCInputError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
