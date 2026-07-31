#!/usr/bin/env python3
"""Extract one publisher-canonical transcript per gene for JCVI.

This helper is intended for releases where the publisher explicitly numbers a
canonical transcript (for example ``.t1``) and the selected protein identifiers
close exactly to the GFF3 transcript identifiers and gene count.  It never
guesses a canonical isoform when that exact closure is absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path


class CanonicalInputError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, str | int]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def parse_attributes(text: str, path: Path, line: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise CanonicalInputError(f"{path.name}:{line}: malformed GFF3 attribute")
        key, value = item.split("=", 1)
        if not key or not value or key in attributes:
            raise CanonicalInputError(f"{path.name}:{line}: invalid/duplicate GFF3 attribute")
        attributes[key] = value
    return attributes


def read_selected_proteins(
    path: Path, canonical_pattern: re.Pattern[str]
) -> tuple[list[tuple[str, str]], set[str]]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    selected: list[tuple[str, str]] = []
    identifiers: set[str] = set()
    current_id: str | None = None
    current_sequence: list[str] = []

    def finish() -> None:
        nonlocal current_id, current_sequence
        if current_id is None:
            return
        sequence = "".join(current_sequence)
        if not sequence:
            raise CanonicalInputError(f"{path.name}: empty protein {current_id!r}")
        if canonical_pattern.search(current_id):
            if current_id in identifiers:
                raise CanonicalInputError(f"{path.name}: duplicate protein ID {current_id!r}")
            identifiers.add(current_id)
            selected.append((current_id, sequence))

    try:
        with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, start=1):
                if raw.startswith(">"):
                    finish()
                    header = raw[1:].strip()
                    current_id = header.split()[0] if header else ""
                    if not current_id:
                        raise CanonicalInputError(
                            f"{path.name}:{line_number}: empty protein identifier"
                        )
                    current_sequence = []
                else:
                    if current_id is None and raw.strip():
                        raise CanonicalInputError(
                            f"{path.name}:{line_number}: sequence before FASTA header"
                        )
                    current_sequence.append("".join(raw.split()))
        finish()
    except (OSError, UnicodeError) as error:
        raise CanonicalInputError(f"Cannot read protein FASTA: {error}") from error
    if not selected:
        raise CanonicalInputError("Canonical protein pattern selected zero records")
    return selected, identifiers


def read_gff_transcripts(
    path: Path, selected_ids: set[str]
) -> tuple[list[tuple[str, int, int, str, str]], int]:
    transcript_rows: dict[str, tuple[str, int, int, str, str]] = {}
    genes: set[str] = set()
    selected_parents: set[str] = set()
    try:
        with path.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw in enumerate(handle, start=1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n\r").split("\t")
                if len(fields) != 9:
                    raise CanonicalInputError(
                        f"{path.name}:{line_number}: GFF3 row must have nine fields"
                    )
                seqid, _, feature, start_text, end_text, _, strand, _, attr_text = fields
                attributes = parse_attributes(attr_text, path, line_number)
                if feature == "gene":
                    gene_id = attributes.get("ID")
                    if not gene_id or gene_id in genes:
                        raise CanonicalInputError(
                            f"{path.name}:{line_number}: missing/duplicate gene ID"
                        )
                    genes.add(gene_id)
                if feature not in {"mRNA", "transcript"}:
                    continue
                transcript_id = attributes.get("ID")
                if transcript_id not in selected_ids:
                    continue
                parent = attributes.get("Parent", "")
                if not parent or "," in parent:
                    raise CanonicalInputError(
                        f"{path.name}:{line_number}: selected transcript lacks one parent gene"
                    )
                if transcript_id in transcript_rows or parent in selected_parents:
                    raise CanonicalInputError(
                        f"{path.name}:{line_number}: selected transcript/gene is not one-to-one"
                    )
                try:
                    start, end = int(start_text), int(end_text)
                except ValueError as error:
                    raise CanonicalInputError(
                        f"{path.name}:{line_number}: non-integer transcript coordinate"
                    ) from error
                if not seqid or start < 1 or end < start or strand not in {"+", "-"}:
                    raise CanonicalInputError(
                        f"{path.name}:{line_number}: invalid selected transcript interval"
                    )
                transcript_rows[transcript_id] = (seqid, start, end, strand, parent)
                selected_parents.add(parent)
    except (OSError, UnicodeError) as error:
        raise CanonicalInputError(f"Cannot read GFF3: {error}") from error
    if set(transcript_rows) != selected_ids:
        missing = sorted(selected_ids.difference(transcript_rows))[:5]
        extra = sorted(set(transcript_rows).difference(selected_ids))[:5]
        raise CanonicalInputError(
            f"Canonical protein/GFF transcript IDs differ: missing={missing}; extra={extra}"
        )
    if len(selected_parents) != len(genes) or selected_parents != genes:
        raise CanonicalInputError(
            "Canonical transcripts do not close one-to-one to the complete GFF gene set"
        )
    rows = [
        (seqid, start, end, transcript_id, strand)
        for transcript_id, (seqid, start, end, strand, _) in transcript_rows.items()
    ]
    rows.sort(key=lambda row: (natural_key(row[0]), row[1], row[3]))
    return rows, len(genes)


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value))


def run(args: argparse.Namespace) -> Path:
    protein = args.protein.expanduser().resolve()
    gff = args.gff.expanduser().resolve()
    for path, label in ((protein, "protein"), (gff, "GFF3")):
        if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
            raise CanonicalInputError(f"{label} input is not a non-empty regular file")
    try:
        pattern = re.compile(args.canonical_id_regex)
    except re.error as error:
        raise CanonicalInputError(f"Invalid canonical ID regex: {error}") from error
    proteins, protein_ids = read_selected_proteins(protein, pattern)
    coordinates, gene_count = read_gff_transcripts(gff, protein_ids)
    if len(proteins) != gene_count:
        raise CanonicalInputError("Canonical protein count does not equal GFF gene count")

    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise CanonicalInputError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        protein_out = staging / f"{args.sample_id}.canonical.protein.faa"
        bed_out = staging / f"{args.sample_id}.canonical.bed"
        coords_out = staging / f"{args.sample_id}.canonical.coords.tsv"
        with protein_out.open("w", encoding="utf-8") as handle:
            for identifier, sequence in proteins:
                handle.write(f">{identifier}\n")
                for start in range(0, len(sequence), 60):
                    handle.write(sequence[start : start + 60] + "\n")
        with bed_out.open("w", encoding="utf-8") as bed, coords_out.open(
            "w", encoding="utf-8", newline=""
        ) as coords:
            writer = csv.writer(coords, delimiter="\t", lineterminator="\n")
            for seqid, start, end, identifier, strand in coordinates:
                bed.write(f"{seqid}\t{start - 1}\t{end}\t{identifier}\t0\t{strand}\n")
                writer.writerow((identifier, seqid, start, end, strand))
        manifest = {
            "schema_version": 1,
            "status": "PASS",
            "sample_id": args.sample_id,
            "selection_rule": {
                "canonical_id_regex": args.canonical_id_regex,
                "require_one_selected_transcript_per_complete_gff_gene": True,
            },
            "inputs": {"protein": binding(protein), "gff": binding(gff)},
            "counts": {
                "gff_genes": gene_count,
                "selected_transcripts": len(coordinates),
                "selected_proteins": len(proteins),
            },
            "outputs": {
                "protein": binding(protein_out),
                "bed": binding(bed_out),
                "coordinates": binding(coords_out),
            },
            "checks": {
                "canonical_protein_ids_unique": True,
                "protein_gff_transcript_id_identity": True,
                "one_selected_transcript_per_gene": True,
                "complete_gene_closure": True,
            },
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.rename(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--protein", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--canonical-id-regex", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    return p


def main() -> int:
    try:
        output = run(parser().parse_args())
        print(f"PASS\t{output}")
        return 0
    except (CanonicalInputError, OSError, UnicodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
