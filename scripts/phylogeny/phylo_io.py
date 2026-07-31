#!/usr/bin/env python3
"""Small, dependency-free helpers for the phylogeny refactor scripts.

These helpers deliberately avoid Biopython so that the validation commands can
run in a minimal Python 3 environment.  They are not a replacement for a
phylogenetic library; they only implement the narrow FASTA/TSV operations used
by the companion scripts.
"""

from __future__ import annotations

import csv
import os
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping


STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class DataError(ValueError):
    """Raised for an input condition that makes an analysis unsafe to continue."""


def normalize_gene_id(value: str) -> str:
    """Return the legacy-compatible FASTA ID used after colon normalisation.

    The old workflow replaced ':' with '_' in CDS headers before matching them
    to protein alignments.  Keeping that normalisation explicit lets the new
    scripts validate the historical files without silently guessing elsewhere.
    """

    return value.strip().split()[0].replace(":", "_")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(id, full_header, sequence)`` records from a FASTA file.

    Empty lines are ignored.  Duplicate IDs are reported by :func:`read_fasta`
    rather than here so callers that need streaming can still use this helper.
    """

    header: str | None = None
    sequence_parts: list[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    record_id = normalize_gene_id(header[1:])
                    yield record_id, header[1:], "".join(sequence_parts).upper()
                header = line
                sequence_parts = []
            elif header is None:
                raise DataError(f"{path}: sequence encountered before a FASTA header")
            else:
                sequence_parts.append(line.replace(" ", ""))
    if header is not None:
        record_id = normalize_gene_id(header[1:])
        yield record_id, header[1:], "".join(sequence_parts).upper()


def read_fasta(path: str | Path) -> "OrderedDict[str, tuple[str, str]]":
    """Read FASTA records in input order and reject duplicate identifiers."""

    records: "OrderedDict[str, tuple[str, str]]" = OrderedDict()
    for record_id, header, sequence in iter_fasta(path):
        if record_id in records:
            raise DataError(f"{path}: duplicate FASTA ID after normalisation: {record_id}")
        if not sequence:
            raise DataError(f"{path}: empty sequence for {record_id}")
        records[record_id] = (header, sequence)
    if not records:
        raise DataError(f"{path}: no FASTA records found")
    return records


def validate_equal_lengths(records: Mapping[str, tuple[str, str]], path: str | Path) -> int:
    """Return alignment length or raise if records are not equally aligned."""

    lengths = {len(sequence) for _, sequence in records.values()}
    if len(lengths) != 1:
        details = ", ".join(str(value) for value in sorted(lengths)[:8])
        raise DataError(f"{path}: unequal alignment lengths ({details})")
    return next(iter(lengths))


def read_tsv(path: str | Path, required_columns: Iterable[str] = ()) -> list[dict[str, str]]:
    """Read a UTF-8 tab-separated table and verify its header."""

    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise DataError(f"{path}: TSV has no header")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise DataError(f"{path}: missing required columns: {', '.join(missing)}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            if not any(cleaned.values()):
                continue
            cleaned["__line__"] = str(line_number)
            rows.append(cleaned)
    return rows


def bool_value(value: str) -> bool:
    """Parse a strict but friendly Boolean TSV value."""

    normalised = value.strip().lower()
    if normalised in {"true", "yes", "1", "y"}:
        return True
    if normalised in {"false", "no", "0", "n", "", "na"}:
        return False
    raise DataError(f"expected true/false, found {value!r}")


@dataclass(frozen=True)
class TerminalManifest:
    """Normalised terminal selection from a taxon registry or terminal table.

    ``selected`` contains only rows explicitly enabled for the current species
    tree.  Candidate and excluded rows remain available for metadata reporting
    but are never silently promoted into an alignment, tree, or CAFE matrix.
    """

    schema: str
    all_rows: list[dict[str, str]]
    selected: list[dict[str, str]]
    candidates: list[dict[str, str]]
    excluded: list[dict[str, str]]


def _first_field(row: Mapping[str, str], names: Iterable[str]) -> tuple[str, str] | None:
    for name in names:
        if name in row and row[name].strip():
            return name, row[name].strip()
    return None


def _tree_selection(value: str, path: str | Path, line: str) -> str:
    normalised = value.strip().lower()
    if normalised in {"true", "yes", "1", "y", "selected", "active", "include"}:
        return "selected"
    if normalised in {
        "candidate",
        "pending",
        "pending_qc",
        "pending_design",
        "parental_lineage_model",
    }:
        return "candidate"
    if normalised in {"false", "no", "0", "n", "", "excluded", "exclude"}:
        return "excluded"
    raise DataError(
        f"{path}:{line}: include_species_tree must be true, false, candidate, or pending; "
        f"found {value!r}"
    )


def load_terminal_manifest(path: str | Path) -> TerminalManifest:
    """Load either ``taxa.tsv`` or a generic phylogeny-terminal manifest.

    The taxon-registry schema uses ``taxon_id``, ``biological_species`` and
    ``role``.  A terminal manifest may instead use ``terminal_id`` or the
    backwards-compatible ``sample_id`` and must group technical terminals via
    ``biological_species``, ``biological_species_id`` or ``grouping``.  Stable
    FASTA stems and tree labels default to the terminal identifier only for the
    compact taxon-registry form.
    """

    raw_rows = read_tsv(path)
    if not raw_rows:
        raise DataError(f"{path}: terminal manifest is empty")
    columns = set(raw_rows[0])
    schema = "taxon_registry" if "taxon_id" in columns else "terminal_manifest"

    normalised_rows: list[dict[str, str]] = []
    for raw in raw_rows:
        line = raw["__line__"]
        identifier = _first_field(raw, ("terminal_id", "sample_id", "taxon_id", "assembly_unit_id"))
        if identifier is None:
            raise DataError(
                f"{path}:{line}: missing terminal identifier; expected terminal_id, sample_id, "
                "taxon_id, or assembly_unit_id"
            )
        _, terminal_id = identifier
        species = _first_field(raw, ("biological_species", "biological_species_id", "grouping"))
        if species is None:
            raise DataError(
                f"{path}:{line}: missing biological-species grouping; expected biological_species, "
                "biological_species_id, or grouping"
            )
        _, biological_species = species
        grouping = _first_field(raw, ("grouping", "biological_species_id", "biological_species"))
        assert grouping is not None

        source = _first_field(raw, ("source_fasta_stem",))
        label = _first_field(raw, ("canonical_tree_label", "tree_label"))
        if schema == "terminal_manifest" and source is None:
            raise DataError(f"{path}:{line}: terminal manifest has no source_fasta_stem")
        if schema == "terminal_manifest" and label is None:
            raise DataError(f"{path}:{line}: terminal manifest has no canonical_tree_label")

        selection = _tree_selection(raw.get("include_species_tree", "true"), path, line)
        root_value = raw.get("is_root_outgroup", "false")
        # Parse now so malformed values fail even when a row is only a candidate.
        is_root = bool_value(root_value)
        row = dict(raw)
        row.update(
            {
                "sample_id": terminal_id,
                "terminal_id": terminal_id,
                "biological_species": biological_species,
                "biological_species_id": grouping[1],
                "grouping": grouping[1],
                "source_fasta_stem": source[1] if source else terminal_id,
                "canonical_tree_label": label[1] if label else terminal_id,
                "is_root_outgroup": "true" if is_root else "false",
                "identity_status": raw.get("identity_status", raw.get("current_asset_status", "")),
                "tree_selection_status": selection,
            }
        )
        normalised_rows.append(row)

    for field in ("sample_id", "source_fasta_stem", "canonical_tree_label"):
        seen: dict[str, int] = {}
        for row in normalised_rows:
            value = row[field]
            if not value:
                raise DataError(f"{path}:{row['__line__']}: empty {field}")
            if value in seen:
                raise DataError(
                    f"{path}:{row['__line__']}: duplicate {field}={value!r}; "
                    f"first used at line {seen[value]}"
                )
            seen[value] = int(row["__line__"])

    return TerminalManifest(
        schema=schema,
        all_rows=normalised_rows,
        selected=[row for row in normalised_rows if row["tree_selection_status"] == "selected"],
        candidates=[row for row in normalised_rows if row["tree_selection_status"] == "candidate"],
        excluded=[row for row in normalised_rows if row["tree_selection_status"] == "excluded"],
    )


def load_samples(path: str | Path) -> list[dict[str, str]]:
    """Load terminals explicitly selected for the current species-tree run."""

    manifest = load_terminal_manifest(path)
    if not manifest.selected:
        raise DataError(f"{path}: no rows have include_species_tree=true")
    return manifest.selected


def fasta_source_stem(value: str) -> str:
    """Remove one conventional FASTA extension from an OrthoFinder label."""

    for suffix in (".fasta", ".faa", ".fas", ".fa"):
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def source_label_to_sample(samples: Iterable[dict[str, str]]) -> dict[str, str]:
    """Map the OrthoFinder source label (without .fa) to stable sample_id."""

    mapping: dict[str, str] = {}
    for row in samples:
        stem = fasta_source_stem(row["source_fasta_stem"])
        if stem in mapping and mapping[stem] != row["sample_id"]:
            raise DataError(f"source label maps to multiple samples: {stem}")
        mapping[stem] = row["sample_id"]
    return mapping


def sample_by_id(samples: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["sample_id"]: row for row in samples}


def load_species_id_map(path: str | Path) -> dict[str, str]:
    """Map OrthoFinder numeric species prefixes to source FASTA labels."""

    mapping: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if ": " not in line:
                raise DataError(f"{path}:{line_number}: malformed SpeciesIDs.txt row")
            numeric_id, source_label = line.split(": ", 1)
            if not numeric_id.isdigit() or not source_label:
                raise DataError(f"{path}:{line_number}: malformed SpeciesIDs.txt row")
            if numeric_id in mapping:
                raise DataError(f"{path}:{line_number}: duplicate numeric species ID {numeric_id}")
            mapping[numeric_id] = source_label
    if not mapping:
        raise DataError(f"{path}: no species mappings found")
    return mapping


def load_sequence_id_map(
    path: str | Path,
    needed_ids: set[str] | None = None,
    species_ids: str | Path | None = None,
) -> dict[str, str]:
    """Map normalised protein IDs to their OrthoFinder source FASTA labels.

    Older OrthoFinder rows embed ``[Species_label]``.  OrthoFinder 3 writes
    ``0_0: gene_id`` and stores source labels in ``SpeciesIDs.txt``.  Both
    layouts are accepted and cross-checked when both label sources exist.
    """

    sequence_path = Path(path)
    species_path = Path(species_ids) if species_ids is not None else sequence_path.with_name("SpeciesIDs.txt")
    numeric_species: dict[str, str] | None = None
    mapping: dict[str, str] = {}
    with sequence_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if ": " not in line:
                raise DataError(f"{sequence_path}:{line_number}: malformed SequenceIDs.txt row")
            internal_id, rest = line.split(": ", 1)
            species_prefix = internal_id.split("_", 1)[0]
            if not species_prefix.isdigit() or not rest:
                raise DataError(f"{sequence_path}:{line_number}: malformed SequenceIDs.txt row")
            if " [" in rest and rest.endswith("]"):
                gene_id, source_label = rest.rsplit(" [", 1)
                source_label = source_label[:-1]
                if species_path.exists():
                    if numeric_species is None:
                        numeric_species = load_species_id_map(species_path)
                    expected_source = numeric_species.get(species_prefix)
                    if expected_source != source_label:
                        raise DataError(
                            f"{sequence_path}:{line_number}: embedded source {source_label!r} "
                            f"does not match {expected_source!r} in {species_path}"
                        )
            else:
                gene_id = rest
                if numeric_species is None:
                    if not species_path.is_file():
                        raise DataError(
                            f"{sequence_path}:{line_number}: OrthoFinder 3 row lacks an embedded "
                            f"source label and {species_path} is missing"
                        )
                    numeric_species = load_species_id_map(species_path)
                source_label = numeric_species.get(species_prefix)
                if source_label is None:
                    raise DataError(
                        f"{sequence_path}:{line_number}: species prefix {species_prefix} "
                        f"is absent from {species_path}"
                    )
            normalised = normalize_gene_id(gene_id)
            if needed_ids is not None and normalised not in needed_ids:
                continue
            old = mapping.get(normalised)
            if old is not None and old != source_label:
                raise DataError(
                    f"{sequence_path}:{line_number}: {normalised} maps to both {old} and {source_label}"
                )
            mapping[normalised] = source_label
    if needed_ids is not None:
        missing = sorted(needed_ids.difference(mapping))
        if missing:
            preview = ", ".join(missing[:10])
            suffix = " ..." if len(missing) > 10 else ""
            raise DataError(
                f"{sequence_path}: {len(missing)} alignment IDs are absent from SequenceIDs.txt: "
                f"{preview}{suffix}"
            )
    return mapping


def resolve_record_samples(
    record_ids: Iterable[str],
    sequence_to_source: Mapping[str, str],
    source_to_sample: Mapping[str, str],
) -> dict[str, str]:
    """Return ``record_id -> sample_id`` and reject unmappable source labels."""

    result: dict[str, str] = {}
    for record_id in record_ids:
        source_label = sequence_to_source.get(record_id)
        if source_label is None:
            raise DataError(f"no source species mapping for alignment ID {record_id}")
        sample_id = source_to_sample.get(fasta_source_stem(source_label))
        if sample_id is None:
            raise DataError(
                f"source label {source_label!r} for {record_id} is absent from samples.tsv"
            )
        result[record_id] = sample_id
    return result


def translate_standard(dna: str) -> str:
    """Translate standard-code DNA; ambiguous codons become ``X``."""

    cleaned = dna.upper().replace("U", "T")
    if len(cleaned) % 3:
        raise DataError(f"CDS length is not divisible by 3 ({len(cleaned)} bp)")
    amino_acids: list[str] = []
    for start in range(0, len(cleaned), 3):
        amino_acids.append(STANDARD_CODE.get(cleaned[start : start + 3], "X"))
    return "".join(amino_acids)


def safe_output_path(path: str | Path, overwrite: bool) -> Path:
    """Reject accidental overwrites and ensure the parent output directory exists."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        raise DataError(f"refusing to overwrite existing output: {output} (use --overwrite)")
    return output


def atomic_write_text(path: str | Path, text: str, overwrite: bool = False) -> None:
    """Write a text result atomically after explicit overwrite checking."""

    output = safe_output_path(path, overwrite)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent, text=True)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary_name, output)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def report_tsv(rows: Iterable[Mapping[str, object]], path: str | Path | None = None) -> str:
    """Render report dictionaries as TSV, optionally persist them atomically."""

    ordered_rows = list(rows)
    columns: list[str] = []
    for row in ordered_rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    if not columns:
        columns = ["check", "severity", "status", "detail"]
    lines = ["\t".join(columns)]
    for row in ordered_rows:
        lines.append(
            "\t".join(str(row.get(column, "")).replace("\t", " ").replace("\n", " ") for column in columns)
        )
    text = "\n".join(lines) + "\n"
    if path is not None and str(path) != "-":
        atomic_write_text(path, text)
    return text
