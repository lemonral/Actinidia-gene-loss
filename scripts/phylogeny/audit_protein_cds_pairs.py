#!/usr/bin/env python3
"""Audit matched protein/CDS FASTA pairs before orthology and codon trees.

The input manifest stores paths relative to a private data root.  The audit is
read-only: it writes a summary and one rejected-ID table per terminal into a
new output directory.  A rejected CDS is not allowed into codon
back-translation, but its protein can still be used for OrthoFinder when the
manifest explicitly permits that use.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

from phylo_io import DataError, read_fasta, read_tsv, translate_standard


REQUIRED_COLUMNS = (
    "terminal_id",
    "protein_path",
    "cds_path",
    "use_for_orthofinder",
    "use_for_codon_tree",
)


@dataclass(frozen=True)
class AuditSummary:
    terminal_id: str
    protein_path: str
    cds_path: str
    protein_sha256: str
    cds_sha256: str
    protein_records: int
    cds_records: int
    shared_ids: int
    protein_only_ids: int
    cds_only_ids: int
    frame_failures: int
    internal_stop_failures: int
    translation_failures: int
    codon_eligible_ids: int
    use_for_orthofinder: str
    use_for_codon_tree: str
    codon_gate: str


def parse_bool(value: str, *, field: str, terminal_id: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "yes", "1"}:
        return True
    if normalized in {"false", "no", "0"}:
        return False
    raise DataError(f"{terminal_id}: {field} must be true or false, found {value!r}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def first_translation_difference(protein: str, translated: str) -> int | None:
    protein = protein.rstrip("*")
    translated = translated.rstrip("*")
    if len(protein) != len(translated):
        return min(len(protein), len(translated)) + 1
    for index, (left, right) in enumerate(zip(protein, translated), start=1):
        if left == right or "X" in {left, right}:
            continue
        return index
    return None


def resolve_input(data_root: Path, value: str, terminal_id: str, field: str) -> Path:
    candidate = Path(value)
    path = candidate if candidate.is_absolute() else data_root / candidate
    if not path.is_file():
        raise DataError(f"{terminal_id}: {field} is not a file: {path}")
    return path


def atomic_write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def audit_pair(row: dict[str, str], data_root: Path) -> tuple[AuditSummary, list[dict[str, object]]]:
    terminal_id = row["terminal_id"]
    protein_path = resolve_input(data_root, row["protein_path"], terminal_id, "protein_path")
    cds_path = resolve_input(data_root, row["cds_path"], terminal_id, "cds_path")
    use_orthofinder = parse_bool(
        row["use_for_orthofinder"], field="use_for_orthofinder", terminal_id=terminal_id
    )
    use_codon = parse_bool(
        row["use_for_codon_tree"], field="use_for_codon_tree", terminal_id=terminal_id
    )

    proteins = read_fasta(protein_path)
    cds = read_fasta(cds_path)
    protein_ids = set(proteins)
    cds_ids = set(cds)
    shared = protein_ids & cds_ids
    rejected: list[dict[str, object]] = []

    for record_id in sorted(protein_ids - cds_ids):
        rejected.append({"terminal_id": terminal_id, "record_id": record_id, "reason": "cds_missing", "detail": ""})
    for record_id in sorted(cds_ids - protein_ids):
        rejected.append({"terminal_id": terminal_id, "record_id": record_id, "reason": "protein_missing", "detail": ""})

    frame_failures = internal_stops = translation_failures = 0
    eligible = 0
    for record_id in sorted(shared):
        protein = proteins[record_id][1]
        nucleotide = cds[record_id][1]
        if len(nucleotide) % 3:
            frame_failures += 1
            rejected.append(
                {
                    "terminal_id": terminal_id,
                    "record_id": record_id,
                    "reason": "cds_length_not_divisible_by_three",
                    "detail": str(len(nucleotide)),
                }
            )
            continue
        translated = translate_standard(nucleotide)
        internal_stop = "*" in translated.rstrip("*")
        difference = first_translation_difference(protein, translated)
        if internal_stop:
            internal_stops += 1
            rejected.append(
                {
                    "terminal_id": terminal_id,
                    "record_id": record_id,
                    "reason": "internal_stop",
                    "detail": "",
                }
            )
        if difference is not None:
            translation_failures += 1
            rejected.append(
                {
                    "terminal_id": terminal_id,
                    "record_id": record_id,
                    "reason": "protein_cds_translation_mismatch",
                    "detail": str(difference),
                }
            )
        if not internal_stop and difference is None:
            eligible += 1

    codon_gate = "PASS" if (not use_codon or eligible > 0) else "FAIL"
    summary = AuditSummary(
        terminal_id=terminal_id,
        protein_path=row["protein_path"],
        cds_path=row["cds_path"],
        protein_sha256=sha256_file(protein_path),
        cds_sha256=sha256_file(cds_path),
        protein_records=len(proteins),
        cds_records=len(cds),
        shared_ids=len(shared),
        protein_only_ids=len(protein_ids - cds_ids),
        cds_only_ids=len(cds_ids - protein_ids),
        frame_failures=frame_failures,
        internal_stop_failures=internal_stops,
        translation_failures=translation_failures,
        codon_eligible_ids=eligible,
        use_for_orthofinder=str(use_orthofinder).lower(),
        use_for_codon_tree=str(use_codon).lower(),
        codon_gate=codon_gate,
    )
    return summary, rejected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    output_dir = Path(args.output_dir)
    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        raise SystemExit(f"refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(manifest_path, REQUIRED_COLUMNS)
    terminal_ids = [row["terminal_id"] for row in rows]
    if len(terminal_ids) != len(set(terminal_ids)):
        raise SystemExit("manifest has duplicate terminal_id values")

    summaries: list[AuditSummary] = []
    output_files: list[Path] = []
    for row in rows:
        try:
            summary, rejected = audit_pair(row, Path(args.data_root))
        except DataError as error:
            raise SystemExit(str(error)) from error
        summaries.append(summary)
        rejected_path = output_dir / f"{summary.terminal_id}.rejected_ids.tsv"
        atomic_write_tsv(
            rejected_path,
            ["terminal_id", "record_id", "reason", "detail"],
            rejected,
        )
        output_files.append(rejected_path)

    summary_path = output_dir / "protein_cds_pair_audit.tsv"
    summary_rows = [asdict(item) for item in summaries]
    atomic_write_tsv(summary_path, list(summary_rows[0]), summary_rows)
    output_files.append(summary_path)

    provenance = {
        "schema_version": 1,
        "manifest_basename": manifest_path.name,
        "manifest_sha256": sha256_file(manifest_path),
        "terminal_count": len(summaries),
        "outputs": [
            {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(output_files)
        ],
    }
    provenance_path = output_dir / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(summary_path)
    return 2 if any(item.codon_gate == "FAIL" for item in summaries) else 0


if __name__ == "__main__":
    raise SystemExit(main())
