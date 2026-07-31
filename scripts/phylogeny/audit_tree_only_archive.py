#!/usr/bin/env python3
"""Audit a published tree-only ZIP before it is used in a new alignment.

The audit is deliberately descriptive.  It reports where an exact taxon label
occurs, whether that label has an assembled target-locus record, and how much
non-missing sequence it contributes to a named concatenated alignment.  It
does not infer that a published terminal is compatible with a newly generated
orthologue matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import BinaryIO, Iterable
from zipfile import ZipFile, ZipInfo


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def member_contains(
    handle: BinaryIO,
    needle: bytes,
    *,
    chunk_size: int = 1024 * 1024,
) -> bool:
    """Search a binary member without assuming it fits in memory."""

    if not needle:
        raise ValueError("needle must not be empty")
    overlap = b""
    needle_lower = needle.lower()
    keep = max(0, len(needle) - 1)
    while chunk := handle.read(chunk_size):
        data = overlap + chunk
        if needle_lower in data.lower():
            return True
        overlap = data[-keep:] if keep else b""
    return False


def fasta_headers(lines: Iterable[str]) -> Iterable[str]:
    for line in lines:
        if line.startswith(">"):
            yield line[1:].strip().split(maxsplit=1)[0]


def fasta_record(lines: Iterable[str], target_id: str) -> str | None:
    found = False
    sequence: list[str] = []
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if found:
                break
            found = line[1:].split(maxsplit=1)[0] == target_id
        elif found:
            sequence.append(line)
    return "".join(sequence) if found else None


def nonmissing_blocks(sequence: str, missing: frozenset[str]) -> list[list[int]]:
    positions = [
        index + 1
        for index, character in enumerate(sequence)
        if character.upper() not in missing
    ]
    if not positions:
        return []
    blocks: list[list[int]] = []
    start = end = positions[0]
    for position in positions[1:]:
        if position == end + 1:
            end = position
        else:
            blocks.append([start, end])
            start = end = position
    blocks.append([start, end])
    return blocks


def audit_archive(
    archive: Path,
    taxon_label: str,
    assembled_prefix: str,
    concat_member: str,
    missing_characters: str = "-?NX.",
) -> dict[str, object]:
    missing = frozenset(character.upper() for character in missing_characters)
    label_bytes = taxon_label.encode("utf-8")
    matching_members: list[str] = []
    assembled_matching_headers: list[dict[str, str]] = []

    with ZipFile(archive) as bundle:
        members = [info for info in bundle.infolist() if not info.is_dir()]
        names = {info.filename for info in members}
        if concat_member not in names:
            raise ValueError(f"concatenated alignment member is absent: {concat_member}")

        for info in members:
            with bundle.open(info) as handle:
                if member_contains(handle, label_bytes):
                    matching_members.append(info.filename)

            if info.filename.startswith(assembled_prefix):
                with bundle.open(info) as handle:
                    text = (
                        line.decode("utf-8", "replace")
                        for line in handle
                    )
                    for header in fasta_headers(text):
                        if header == taxon_label:
                            assembled_matching_headers.append(
                                {"member": info.filename, "header": header}
                            )

        with bundle.open(concat_member) as handle:
            text = (line.decode("utf-8", "replace") for line in handle)
            sequence = fasta_record(text, taxon_label)

    blocks = nonmissing_blocks(sequence or "", missing)
    nonmissing_sites = sum(end - start + 1 for start, end in blocks)
    concat_stats = {
        "member": concat_member,
        "record_found": sequence is not None,
        "alignment_length": len(sequence or ""),
        "nonmissing_sites": nonmissing_sites,
        "nonmissing_fraction": (
            nonmissing_sites / len(sequence) if sequence else None
        ),
        "first_nonmissing_1based": blocks[0][0] if blocks else None,
        "last_nonmissing_1based": blocks[-1][1] if blocks else None,
        "nonmissing_blocks_1based_inclusive": blocks,
        "missing_characters": "".join(sorted(missing)),
    }

    return {
        "schema_version": "1.0.0",
        "archive": str(archive),
        "archive_bytes": archive.stat().st_size,
        "archive_sha256": sha256_file(archive),
        "zip_member_count": len(members),
        "taxon_label": taxon_label,
        "matching_members": sorted(matching_members),
        "assembled_prefix": assembled_prefix,
        "assembled_exact_header_matches": assembled_matching_headers,
        "assembled_exact_header_match_count": len(assembled_matching_headers),
        "concatenated_alignment": concat_stats,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--taxon-label", required=True)
    parser.add_argument("--assembled-prefix", required=True)
    parser.add_argument("--concat-member", required=True)
    parser.add_argument("--missing-characters", default="-?NX.")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = audit_archive(
        archive=args.archive.resolve(),
        taxon_label=args.taxon_label,
        assembled_prefix=args.assembled_prefix,
        concat_member=args.concat_member,
        missing_characters=args.missing_characters,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
