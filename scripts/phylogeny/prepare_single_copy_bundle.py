#!/usr/bin/env python3
"""Build a fail-closed protein/CDS bundle from OrthoFinder single-copy loci.

This workflow validates every original single-copy protein through
SpeciesIDs/SequenceIDs, binds it to the audited terminal-specific CDS file,
and publishes only after all loci pass together.  OrthoFinder's internal gene-
tree alignments are intentionally not reused because its reconciliation stage
may retain only a homologous subsequence of a long protein; the published
bundle is the input for a fresh, full-length MAFFT run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from phylo_io import (
    DataError,
    fasta_source_stem,
    iter_fasta,
    load_samples,
    load_sequence_id_map,
    read_fasta,
    read_tsv,
    source_label_to_sample,
    translate_standard,
)


PAIR_COLUMNS = (
    "terminal_id",
    "cds_path",
    "use_for_codon_tree",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_true(value: str, context: str) -> bool:
    normalised = value.strip().lower()
    if normalised in {"true", "yes", "1"}:
        return True
    if normalised in {"false", "no", "0"}:
        return False
    raise DataError(f"{context}: expected true/false, found {value!r}")


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n{sequence}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--orthofinder-results", required=True, type=Path)
    parser.add_argument("--sequence-pairs", required=True, type=Path)
    parser.add_argument("--terminals", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-loci", type=int, default=479)
    args = parser.parse_args()

    temporary: Path | None = None
    try:
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite existing output directory: {args.output_dir}")
        results = args.orthofinder_results.resolve()
        sequence_ids = results / "WorkingDirectory" / "SequenceIDs.txt"
        species_ids = results / "WorkingDirectory" / "SpeciesIDs.txt"
        sco_list = results / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
        unaligned_dir = results / "Single_Copy_Orthologue_Sequences"
        for required in (sequence_ids, species_ids, sco_list, unaligned_dir):
            if not required.exists():
                raise DataError(f"required OrthoFinder output is missing: {required}")

        loci = [line.strip() for line in sco_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(loci) != args.expected_loci or len(set(loci)) != len(loci):
            raise DataError(
                f"single-copy list has {len(loci)} rows/{len(set(loci))} unique; "
                f"expected {args.expected_loci}"
            )

        samples = load_samples(args.terminals)
        terminal_ids = {row["sample_id"] for row in samples}
        source_to_terminal = source_label_to_sample(samples)
        if len(samples) != len(source_to_terminal):
            raise DataError("terminal manifest does not have one unique source stem per terminal")
        pair_rows = read_tsv(args.sequence_pairs, PAIR_COLUMNS)
        cds_paths: dict[str, Path] = {}
        for row in pair_rows:
            if not parse_true(row["use_for_codon_tree"], f"{args.sequence_pairs}:{row['__line__']}"):
                continue
            terminal_id = row["terminal_id"]
            if terminal_id in cds_paths:
                raise DataError(f"duplicate CDS row for {terminal_id}")
            cds_paths[terminal_id] = (args.data_root / row["cds_path"]).resolve()
        if set(cds_paths) != terminal_ids:
            raise DataError(
                "terminal/CDS set mismatch: missing="
                + ",".join(sorted(terminal_ids - set(cds_paths)))
                + "; extra="
                + ",".join(sorted(set(cds_paths) - terminal_ids))
            )
        for terminal_id, path in cds_paths.items():
            if not path.is_file():
                raise DataError(f"{terminal_id}: CDS FASTA is missing: {path}")

        required_ids: set[str] = set()
        unaligned_by_locus: dict[str, dict[str, tuple[str, str]]] = {}
        for locus in loci:
            path = unaligned_dir / f"{locus}.fa"
            records = read_fasta(path)
            if len(records) != len(samples):
                raise DataError(f"{path}: found {len(records)} records; expected {len(samples)}")
            duplicate = required_ids.intersection(records)
            if duplicate:
                raise DataError(f"protein IDs recur across single-copy loci: {sorted(duplicate)[:5]}")
            required_ids.update(records)
            unaligned_by_locus[locus] = records

        sequence_to_source = load_sequence_id_map(sequence_ids, required_ids, species_ids)
        id_to_terminal: dict[str, str] = {}
        per_terminal_ids = {terminal_id: set() for terminal_id in terminal_ids}
        for record_id, source_label in sequence_to_source.items():
            source_stem = fasta_source_stem(source_label)
            terminal_id = source_to_terminal.get(source_stem)
            if terminal_id is None:
                raise DataError(f"{record_id}: unknown OrthoFinder source {source_label!r}")
            id_to_terminal[record_id] = terminal_id
            per_terminal_ids[terminal_id].add(record_id)

        for locus, records in unaligned_by_locus.items():
            seen_terminals = {id_to_terminal[record_id] for record_id in records}
            if seen_terminals != terminal_ids:
                raise DataError(f"{locus}: terminal closure failed")

        cds_records: dict[str, tuple[str, str]] = {}
        cds_source_sha: dict[str, str] = {}
        for terminal_id in sorted(terminal_ids):
            path = cds_paths[terminal_id]
            needed = per_terminal_ids[terminal_id]
            found: set[str] = set()
            for record_id, _, sequence in iter_fasta(path):
                if record_id not in needed:
                    continue
                if record_id in found or record_id in cds_records:
                    raise DataError(f"{path}: duplicate required CDS ID {record_id}")
                found.add(record_id)
                cds_records[record_id] = (terminal_id, sequence)
            if found != needed:
                missing = sorted(needed - found)
                raise DataError(f"{path}: missing {len(missing)} SCO CDS records: {missing[:5]}")
            cds_source_sha[terminal_id] = sha256_file(path)

        for locus in loci:
            for record_id, (_, protein) in unaligned_by_locus[locus].items():
                terminal_id, cds = cds_records[record_id]
                if len(cds) % 3:
                    raise DataError(f"{terminal_id}/{record_id}: CDS length is not divisible by three")
                if translate_standard(cds).rstrip("*") != protein.rstrip("*"):
                    raise DataError(f"{terminal_id}/{record_id}: protein/CDS translation mismatch")

        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
        protein_output = temporary / "protein_unaligned"
        protein_output.mkdir()
        locus_rows: list[dict[str, object]] = []
        output_checksums: list[tuple[str, int, str]] = []
        for locus in loci:
            destination = protein_output / f"{locus}.fa"
            records = [
                (record_id, sequence.rstrip("*"))
                for record_id, (_, sequence) in unaligned_by_locus[locus].items()
            ]
            write_fasta(destination, records)
            checksum = sha256_file(destination)
            output_checksums.append((str(destination.relative_to(temporary)), destination.stat().st_size, checksum))
            locus_rows.append(
                {
                    "orthogroup": locus,
                    "records": len(records),
                    "protein_unaligned_sha256": checksum,
                }
            )

        combined_cds = temporary / "combined_sco_cds.fa"
        write_fasta(combined_cds, [(record_id, cds_records[record_id][1]) for record_id in sorted(required_ids)])
        output_checksums.append(
            (str(combined_cds.relative_to(temporary)), combined_cds.stat().st_size, sha256_file(combined_cds))
        )
        locus_manifest = temporary / "locus_manifest.tsv"
        with locus_manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(locus_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(locus_rows)
        output_checksums.append(
            (str(locus_manifest.relative_to(temporary)), locus_manifest.stat().st_size, sha256_file(locus_manifest))
        )

        checksum_manifest = temporary / "checksums.tsv"
        with checksum_manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("relative_path\tbytes\tsha256\n")
            for relative_path, size, checksum in sorted(output_checksums):
                handle.write(f"{relative_path}\t{size}\t{checksum}\n")

        binding = {
            "schema_version": 1,
            "status": "PASS",
            "workflow": "prepare_single_copy_bundle",
            "terminals": len(samples),
            "loci": len(loci),
            "protein_cds_pairs": len(required_ids),
            "orthofinder_inputs": {
                "single_copy_list_sha256": sha256_file(sco_list),
                "sequence_ids_sha256": sha256_file(sequence_ids),
                "species_ids_sha256": sha256_file(species_ids),
            },
            "configuration": {
                "sequence_pairs_sha256": sha256_file(args.sequence_pairs),
                "terminals_sha256": sha256_file(args.terminals),
            },
            "cds_source_sha256": cds_source_sha,
            "checksum_manifest_sha256": sha256_file(checksum_manifest),
        }
        (temporary / "bundle_validation.json").write_text(
            json.dumps(binding, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output_dir)
        temporary = None
        print(json.dumps(binding, indent=2, sort_keys=True))
        return 0
    except (DataError, OSError, csv.Error, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
