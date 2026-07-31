#!/usr/bin/env python3
"""Append a declared suffix to FASTA IDs with exact target-ID closure.

This utility is intended for legacy protein/CDS pairs whose sequences match
but whose first-token identifiers differ by one deterministic suffix.  It
never changes sequence bytes and refuses to publish unless the remapped source
ID set equals the expected FASTA ID set exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path


class RemapError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].strip().split(maxsplit=1)[0]
            if not identifier:
                raise RemapError(f"{path}:{line_number}: empty FASTA ID")
            if identifier in seen:
                raise RemapError(f"{path}:{line_number}: duplicate FASTA ID {identifier!r}")
            seen.add(identifier)
            identifiers.append(identifier)
    if not identifiers:
        raise RemapError(f"{path}: no FASTA records")
    return identifiers


def remap(
    source: Path,
    expected: Path,
    output: Path,
    suffix: str,
    *,
    allow_source_extra: bool = False,
) -> dict[str, object]:
    if not suffix or any(character.isspace() for character in suffix):
        raise RemapError("suffix must be non-empty and contain no whitespace")
    if output.exists():
        raise RemapError(f"refusing to overwrite existing output: {output}")

    source_ids = fasta_ids(source)
    expected_ids = fasta_ids(expected)
    remapped_ids = [identifier + suffix for identifier in source_ids]
    if len(remapped_ids) != len(set(remapped_ids)):
        raise RemapError("remapping produced duplicate FASTA IDs")
    missing_ids = sorted(set(expected_ids) - set(remapped_ids))
    extra_ids = sorted(set(remapped_ids) - set(expected_ids))
    if missing_ids or (extra_ids and not allow_source_extra):
        raise RemapError(
            "remapped source IDs do not equal expected IDs; "
            f"missing={missing_ids[:5]!r}; extra={extra_ids[:5]!r}"
        )
    expected_set = set(expected_ids)

    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", dir=output.parent
    )
    try:
        with source.open("r", encoding="utf-8") as source_handle, os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as output_handle:
            record_count = 0
            current_included = False
            for line in source_handle:
                if line.startswith(">"):
                    header = line[1:].rstrip("\r\n")
                    parts = header.split(maxsplit=1)
                    replacement = parts[0] + suffix
                    description = f" {parts[1]}" if len(parts) == 2 else ""
                    current_included = replacement in expected_set
                    if current_included:
                        output_handle.write(f">{replacement}{description}\n")
                        record_count += 1
                else:
                    if current_included:
                        output_handle.write(line)
        os.replace(temporary_name, output)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise

    if set(fasta_ids(output)) != expected_set or record_count != len(expected_ids):
        raise RemapError("published FASTA IDs do not equal the expected ID set")
    extra_digest = hashlib.sha256(
        ("\n".join(extra_ids) + ("\n" if extra_ids else "")).encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": 1,
        "transform": "append_suffix_to_first_token_id",
        "suffix": suffix,
        "record_count": record_count,
        "source_extra_id_count": len(extra_ids),
        "source_extra_ids_sha256": extra_digest,
        "source_sha256": sha256_file(source),
        "expected_ids_fasta_sha256": sha256_file(expected),
        "output_sha256": sha256_file(output),
        "id_set_closure": "PASS",
        "sequence_lines_copied_without_transformation": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--expected-ids-fasta", required=True, type=Path)
    parser.add_argument("--suffix", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument(
        "--allow-source-extra",
        action="store_true",
        help="Filter remapped source IDs absent from the expected FASTA; missing expected IDs remain fatal.",
    )
    args = parser.parse_args()

    if args.provenance.exists():
        raise SystemExit(f"refusing to overwrite existing provenance: {args.provenance}")
    try:
        provenance = remap(
            args.source,
            args.expected_ids_fasta,
            args.output,
            args.suffix,
            allow_source_extra=args.allow_source_extra,
        )
    except RemapError as error:
        raise SystemExit(str(error)) from error
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
