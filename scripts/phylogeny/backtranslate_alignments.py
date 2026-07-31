#!/usr/bin/env python3
"""Create codon alignments from protein alignments with explicit frame checks.

This is a safe replacement for the old ``alignment_protein_to_dna.py`` helper.
It scans a combined CDS FASTA only once, retains just the IDs present in the
selected protein alignments, and refuses to write any alignment if a required
CDS is absent, out of frame, or translates inconsistently.

It is intentionally a *new-run* helper.  It does not assert that it can
recreate the historical alignment byte-for-byte, because the old run did not
preserve a complete command/version/input snapshot.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_io import (
    DataError,
    atomic_write_text,
    iter_fasta,
    read_fasta,
    report_tsv,
    safe_output_path,
    translate_standard,
    validate_equal_lengths,
)


def first_difference(left: str, right: str, allow_x: bool) -> int | None:
    """Return first non-matching amino-acid index, optionally treating X as wildcard."""

    if len(left) != len(right):
        return min(len(left), len(right))
    for index, (a, b) in enumerate(zip(left, right), start=1):
        if a == b or (allow_x and "X" in {a, b}):
            continue
        return index
    return None


def build_codon_alignment(
    protein_records: dict[str, tuple[str, str]],
    cds_records: dict[str, tuple[str, str]],
    allow_x: bool,
) -> tuple[str, list[dict[str, object]]]:
    """Backtranslate one protein alignment after validating every sequence."""

    output: list[tuple[str, str]] = []
    report: list[dict[str, object]] = []
    for record_id, (header, protein_gapped) in protein_records.items():
        if record_id not in cds_records:
            raise DataError(f"{record_id}: CDS is missing")
        _, cds = cds_records[record_id]
        if len(cds) % 3:
            raise DataError(f"{record_id}: CDS length {len(cds)} is not divisible by 3")
        protein = protein_gapped.replace("-", "")
        translated = translate_standard(cds)
        # Annotation FASTAs often omit a terminal stop codon, so one terminal
        # stop is allowed only after all amino acids in the protein were used.
        protein_compare = protein.rstrip("*")
        translated_compare = translated.rstrip("*")
        difference = first_difference(protein_compare, translated_compare, allow_x)
        if difference is not None:
            raise DataError(
                f"{record_id}: protein/CDS translation mismatch at amino acid {difference}; "
                f"protein length={len(protein_compare)}, CDS translation length={len(translated_compare)}"
            )
        nucleotides: list[str] = []
        cursor = 0
        for amino_acid in protein_gapped:
            if amino_acid == "-":
                nucleotides.append("---")
                continue
            codon = cds[cursor : cursor + 3]
            if len(codon) != 3:
                raise DataError(f"{record_id}: CDS ended before protein alignment")
            nucleotides.append(codon)
            cursor += 3
        unused = cds[cursor:]
        if unused and not (len(unused) == 3 and translate_standard(unused) == "*"):
            raise DataError(
                f"{record_id}: {len(unused)} unconsumed CDS bases remain after back-translation"
            )
        output.append((header, "".join(nucleotides)))
        report.append(
            {
                "record_id": record_id,
                "protein_aa_ungapped": len(protein),
                "cds_bp": len(cds),
                "terminal_stop_ignored": "yes" if unused else "no",
            }
        )
    lengths = {len(sequence) for _, sequence in output}
    if len(lengths) != 1:
        raise DataError("back-translated records have unequal alignment lengths")
    text = "".join(f">{header}\n{sequence}\n" for header, sequence in output)
    return text, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protein-dir", required=True, help="directory containing one aligned protein FASTA per orthogroup")
    parser.add_argument("--glob", default="*.fasta", help="glob within --protein-dir (default: *.fasta)")
    parser.add_argument("--cds-fasta", required=True, help="combined CDS FASTA; read once and never modified")
    parser.add_argument("--output-dir", required=True, help="new directory for codon alignments")
    parser.add_argument("--suffix", default="_codon.fa", help="suffix for each output (default: _codon.fa)")
    parser.add_argument("--allow-x", action="store_true", help="treat X as a wildcard during translation comparison")
    parser.add_argument("--overwrite", action="store_true", help="allow replacing outputs created by an earlier run")
    parser.add_argument("--report", default=None, help="optional TSV report path; default is stdout")
    args = parser.parse_args()

    report: list[dict[str, object]] = []
    try:
        protein_paths = sorted(Path(args.protein_dir).glob(args.glob))
        if not protein_paths:
            raise DataError(f"{args.protein_dir}: no files match {args.glob!r}")
        protein_alignments: dict[Path, dict[str, tuple[str, str]]] = {}
        required_ids: set[str] = set()
        for path in protein_paths:
            records = read_fasta(path)
            validate_equal_lengths(records, path)
            protein_alignments[path] = records
            required_ids.update(records)

        # Preflight output names before scanning a potentially very large CDS file.
        output_paths = {
            path: Path(args.output_dir) / f"{path.stem}{args.suffix}" for path in protein_paths
        }
        for output_path in output_paths.values():
            safe_output_path(output_path, args.overwrite)

        cds_records: dict[str, tuple[str, str]] = {}
        for record_id, header, sequence in iter_fasta(args.cds_fasta):
            if record_id not in required_ids:
                continue
            if record_id in cds_records:
                raise DataError(f"{args.cds_fasta}: duplicate needed CDS ID {record_id}")
            cds_records[record_id] = (header, sequence)
        missing = sorted(required_ids.difference(cds_records))
        if missing:
            preview = ", ".join(missing[:10])
            suffix = " ..." if len(missing) > 10 else ""
            raise DataError(f"{args.cds_fasta}: missing {len(missing)} required CDS IDs: {preview}{suffix}")

        rendered: dict[Path, str] = {}
        for path, records in protein_alignments.items():
            text, per_record = build_codon_alignment(records, cds_records, args.allow_x)
            rendered[output_paths[path]] = text
            report.append(
                {
                    "check": "backtranslate_group",
                    "severity": "INFO",
                    "status": "PASS",
                    "detail": f"{len(records)} records passed translation/frame checks",
                    "group": path.name,
                    "output": str(output_paths[path]),
                }
            )
            for record_report in per_record:
                report.append(
                    {
                        "check": "backtranslate_record",
                        "severity": "INFO",
                        "status": "PASS",
                        "detail": "translation/frame consistent",
                        "group": path.name,
                        **record_report,
                    }
                )
        # Validation completed for every group before the first output is written.
        for output_path, text in rendered.items():
            atomic_write_text(output_path, text, overwrite=args.overwrite)
        report.insert(
            0,
            {
                "check": "backtranslate_summary",
                "severity": "INFO",
                "status": "PASS",
                "detail": f"wrote {len(rendered)} codon alignments after global preflight",
            },
        )
        print(report_tsv(report, args.report), end="")
        return 0
    except DataError as error:
        report.append({"check": "backtranslate", "severity": "ERROR", "status": "FAIL", "detail": str(error)})
        print(report_tsv(report, args.report), end="")
        return 2


if __name__ == "__main__":
    sys.exit(main())
