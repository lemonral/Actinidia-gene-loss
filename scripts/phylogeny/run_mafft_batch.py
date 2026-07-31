#!/usr/bin/env python3
"""Run and validate a bounded, atomic MAFFT batch for SCO proteins."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from phylo_io import DataError, read_fasta, validate_equal_lengths


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mafft", required=True, type=Path)
    parser.add_argument("--jobs", required=True, type=int)
    parser.add_argument("--expected-loci", required=True, type=int)
    parser.add_argument("--expected-records", required=True, type=int)
    parser.add_argument("--expected-version", default="7.526")
    args = parser.parse_args()

    temporary: Path | None = None
    try:
        if not 1 <= args.jobs <= 15:
            raise DataError("--jobs must be between 1 and 15")
        mafft = args.mafft.resolve()
        if not mafft.is_file() or not os.access(mafft, os.X_OK):
            raise DataError(f"MAFFT is missing or not executable: {mafft}")
        version_run = subprocess.run(
            [str(mafft), "--version"], text=True, capture_output=True, check=False
        )
        version = (version_run.stdout + version_run.stderr).strip()
        if version_run.returncode != 0 or args.expected_version not in version:
            raise DataError(
                f"MAFFT version probe failed exact expectation {args.expected_version!r}: {version!r}"
            )
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite existing output directory: {args.output_dir}")
        inputs = sorted(args.input_dir.glob("*.fa"))
        if len(inputs) != args.expected_loci:
            raise DataError(f"found {len(inputs)} input loci; expected {args.expected_loci}")
        input_sha: dict[Path, str] = {}
        input_records: dict[Path, dict[str, tuple[str, str]]] = {}
        for path in inputs:
            records = read_fasta(path)
            if len(records) != args.expected_records:
                raise DataError(
                    f"{path}: found {len(records)} records; expected {args.expected_records}"
                )
            for record_id, (_, sequence) in records.items():
                if "-" in sequence or sequence.endswith("*"):
                    raise DataError(f"{path}: non-aligner-ready sequence for {record_id}")
            input_sha[path] = sha256_file(path)
            input_records[path] = records

        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
        alignment_dir = temporary / "protein_alignments"
        log_dir = temporary / "logs"
        alignment_dir.mkdir()
        log_dir.mkdir()

        def run_one(path: Path) -> dict[str, object]:
            output = alignment_dir / path.name
            completed = subprocess.run(
                [str(mafft), "--auto", "--thread", "1", str(path)],
                text=True,
                capture_output=True,
                check=False,
            )
            (log_dir / f"{path.stem}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
            if completed.returncode != 0:
                raise DataError(f"{path.name}: MAFFT exited {completed.returncode}")
            output.write_text(completed.stdout, encoding="utf-8")
            if sha256_file(path) != input_sha[path]:
                raise DataError(f"{path}: input changed during MAFFT")
            aligned = read_fasta(output)
            validate_equal_lengths(aligned, output)
            original = input_records[path]
            if set(aligned) != set(original):
                raise DataError(f"{path.name}: output identifiers changed")
            for record_id, (_, sequence) in aligned.items():
                if sequence.replace("-", "") != original[record_id][1]:
                    raise DataError(f"{path.name}: MAFFT changed residues for {record_id}")
            return {
                "locus": path.stem,
                "records": len(aligned),
                "aligned_aa": validate_equal_lengths(aligned, output),
                "input_sha256": input_sha[path],
                "output_sha256": sha256_file(output),
            }

        rows: list[dict[str, object]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(run_one, path): path for path in inputs}
            for future in concurrent.futures.as_completed(futures):
                rows.append(future.result())
        rows.sort(key=lambda row: str(row["locus"]))
        if len(rows) != args.expected_loci:
            raise DataError("internal error: incomplete MAFFT result set")

        manifest = temporary / "alignment_manifest.tsv"
        columns = ("locus", "records", "aligned_aa", "input_sha256", "output_sha256")
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in rows:
                handle.write("\t".join(str(row[column]) for column in columns) + "\n")
        binding = {
            "schema_version": 1,
            "status": "PASS",
            "workflow": "bounded_mafft_sco_batch",
            "jobs": args.jobs,
            "loci": len(rows),
            "records_per_locus": args.expected_records,
            "mafft_version": version,
            "mafft_sha256": sha256_file(mafft),
            "alignment_manifest_sha256": sha256_file(manifest),
        }
        (temporary / "batch_validation.json").write_text(
            json.dumps(binding, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if sha256_file(mafft) != binding["mafft_sha256"]:
            raise DataError("MAFFT executable changed during the batch")
        os.replace(temporary, args.output_dir)
        temporary = None
        print(json.dumps(binding, indent=2, sort_keys=True))
        return 0
    except (DataError, OSError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
