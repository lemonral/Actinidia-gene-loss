#!/usr/bin/env python3
"""Relabel one-copy alignments with canonical terminal labels, atomically."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

from phylo_io import (
    DataError,
    load_samples,
    load_sequence_id_map,
    read_fasta,
    resolve_record_samples,
    source_label_to_sample,
    validate_equal_lengths,
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument("--glob", default="*.fa")
    parser.add_argument("--sequence-ids", required=True, type=Path)
    parser.add_argument("--species-ids", required=True, type=Path)
    parser.add_argument("--terminals", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-loci", required=True, type=int)
    args = parser.parse_args()

    temporary: Path | None = None
    try:
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite existing output directory: {args.output_dir}")
        paths = sorted(args.alignment_dir.glob(args.glob))
        if len(paths) != args.expected_loci:
            raise DataError(f"found {len(paths)} alignments; expected {args.expected_loci}")
        samples = load_samples(args.terminals)
        sample_order = [row["sample_id"] for row in samples]
        labels = {row["sample_id"]: row["canonical_tree_label"] for row in samples}
        if len(set(labels.values())) != len(labels):
            raise DataError("canonical terminal labels are not unique")
        source_to_sample = source_label_to_sample(samples)
        alignments: dict[Path, dict[str, tuple[str, str]]] = {}
        all_ids: set[str] = set()
        for path in paths:
            records = read_fasta(path)
            validate_equal_lengths(records, path)
            alignments[path] = records
            all_ids.update(records)
        sequence_to_source = load_sequence_id_map(
            args.sequence_ids, all_ids, args.species_ids
        )

        rendered: dict[Path, list[tuple[str, str]]] = {}
        for path, records in alignments.items():
            record_to_sample = resolve_record_samples(records, sequence_to_source, source_to_sample)
            counts = Counter(record_to_sample.values())
            if set(counts) != set(sample_order) or any(value != 1 for value in counts.values()):
                raise DataError(f"{path}: does not contain exactly one sequence per terminal")
            sample_to_record = {sample: record for record, sample in record_to_sample.items()}
            rendered[path] = [
                (labels[sample], records[sample_to_record[sample]][1]) for sample in sample_order
            ]

        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
        rows: list[tuple[str, int, int, str]] = []
        for source in paths:
            destination = temporary / source.name
            records = rendered[source]
            with destination.open("w", encoding="utf-8", newline="\n") as handle:
                for label, sequence in records:
                    handle.write(f">{label}\n{sequence}\n")
            rows.append(
                (source.stem, len(records), len(records[0][1]), sha256_file(destination))
            )
        manifest = temporary / "alignment_manifest.tsv"
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("locus\trecords\taligned_columns\tsha256\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")
        binding = {
            "schema_version": 1,
            "status": "PASS",
            "workflow": "canonical_terminal_alignment_relabelling",
            "loci": len(rows),
            "terminals": len(samples),
            "sequence_ids_sha256": sha256_file(args.sequence_ids),
            "species_ids_sha256": sha256_file(args.species_ids),
            "terminals_sha256": sha256_file(args.terminals),
            "alignment_manifest_sha256": sha256_file(manifest),
        }
        (temporary / "bundle_validation.json").write_text(
            json.dumps(binding, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output_dir)
        temporary = None
        print(json.dumps(binding, indent=2, sort_keys=True))
        return 0
    except (DataError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
