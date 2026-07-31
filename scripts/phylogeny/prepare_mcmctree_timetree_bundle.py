#!/usr/bin/env python3
"""Prepare a checksum-bound MCMCTree bundle from a frozen topology and TimeTree bounds."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Phylo


class PreparationError(RuntimeError):
    pass


IUPAC_NUCLEOTIDES = frozenset("ACGTRYSWKMBDHVN?-.")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise PreparationError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_fasta(path: Path) -> dict[str, str]:
    source = regular(path)
    records: dict[str, list[str]] = {}
    current: str | None = None
    with source.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                if not identifier or identifier in records:
                    raise PreparationError(f"missing or duplicate FASTA id at line {line_number}")
                if any(character.isspace() for character in identifier):
                    raise PreparationError(f"whitespace in FASTA id: {identifier}")
                records[identifier] = []
                current = identifier
            elif current is None:
                raise PreparationError(f"sequence before FASTA header at line {line_number}")
            else:
                sequence = line.upper()
                invalid = set(sequence) - IUPAC_NUCLEOTIDES
                if invalid:
                    raise PreparationError(
                        f"unsupported nucleotide symbols for {current}: {''.join(sorted(invalid))}"
                    )
                records[current].append(sequence)
    joined = {identifier: "".join(chunks) for identifier, chunks in records.items()}
    if not joined or any(not sequence for sequence in joined.values()):
        raise PreparationError("empty FASTA record set or sequence")
    lengths = {len(sequence) for sequence in joined.values()}
    if len(lengths) != 1:
        raise PreparationError("FASTA sequences are not aligned to equal length")
    alignment_length = next(iter(lengths))
    if alignment_length % 3:
        raise PreparationError("codon supermatrix length is not divisible by three")
    return joined


def read_tsv(path: Path) -> list[dict[str, str]]:
    source = regular(path)
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise PreparationError(f"missing TSV header: {source}")
        return list(reader)


def active_status(value: str) -> bool:
    normalized = value.strip().lower()
    return normalized == "pass" or normalized == "active" or normalized.startswith("active_")


def render_calibrated_tree(tree: object, calibration_by_clade: dict[int, tuple[float, float]]) -> str:
    def render(clade: object) -> str:
        if clade.clades:
            body = "(" + ",".join(render(child) for child in clade.clades) + ")"
        else:
            if not clade.name or any(character.isspace() for character in clade.name):
                raise PreparationError("tree tip is missing or contains whitespace")
            body = clade.name
        interval = calibration_by_clade.get(id(clade))
        if interval is not None:
            minimum, maximum = interval
            body += f" 'B({minimum:.10g},{maximum:.10g})'"
        return body

    return render(tree.root) + ";"


def write_partitioned_phylip(path: Path, records: dict[str, str], tip_order: list[str]) -> list[int]:
    if set(tip_order) != set(records) or len(tip_order) != len(records):
        raise PreparationError("tree/FASTA tip closure failed before PHYLIP writing")
    partition_lengths: list[int] = []
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for offset in range(3):
            partition = {identifier: records[identifier][offset::3] for identifier in tip_order}
            lengths = {len(sequence) for sequence in partition.values()}
            if len(lengths) != 1:
                raise PreparationError("codon-position partition lengths disagree")
            length = next(iter(lengths))
            partition_lengths.append(length)
            handle.write(f" {len(tip_order)} {length}\n")
            for identifier in tip_order:
                handle.write(f"{identifier}  {partition[identifier]}\n")
            handle.write("\n")
    return partition_lengths


def control_text(
    *, seed: int, usedata: int, outfile: str, mcmcfile: str, burnin: int, sampfreq: int,
    nsample: int, root_age_upper: float
) -> str:
    if usedata not in {0, 2, 3}:
        raise PreparationError("unsupported MCMCTree usedata mode")
    if seed <= 0 or burnin < 0 or sampfreq <= 0 or nsample <= 0 or root_age_upper <= 0:
        raise PreparationError("invalid MCMCTree control parameter")
    usedata_value = f"{usedata} in.BV" if usedata == 2 else str(usedata)
    return f"""          seed = {seed}
       seqfile = codon_positions.phy
      treefile = calibrated_topology.trees
      mcmcfile = {mcmcfile}
       outfile = {outfile}

         ndata = 3
       seqtype = 0
       usedata = {usedata_value}
         clock = 2
       RootAge = '<{root_age_upper:.10g}'

         model = 4
         alpha = 0.5
         ncatG = 5
     cleandata = 0

       BDparas = 1 1 0.1 multiplicative
   kappa_gamma = 6 2
   alpha_gamma = 1 1
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1

      finetune = 1: .1 .1 .1 .1 .1 .1
         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}
    checkpoint = 1 0.01 mcmctree.ckpt
"""


def probe_mcmctree(path: Path) -> tuple[str, str]:
    executable = regular(path)
    if not os.access(executable, os.X_OK):
        raise PreparationError(f"MCMCTree is not executable: {executable}")
    with tempfile.TemporaryDirectory(prefix="mcmctree-version-") as temporary:
        completed = subprocess.run(
            [str(executable)],
            cwd=temporary,
            input="\n",
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )
    banner_text = completed.stdout + "\n" + completed.stderr
    match = re.search(r"MCMCTREE in paml version ([0-9.]+)", banner_text)
    if not match:
        raise PreparationError("could not recover the MCMCTree version banner")
    return match.group(1), sha256(executable)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--topology-freeze", required=True, type=Path)
    parser.add_argument("--dating-gate", required=True, type=Path)
    parser.add_argument("--constraints", required=True, type=Path)
    parser.add_argument("--mcmctree", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--time-unit-ma", type=float, default=100.0)
    args = parser.parse_args()
    try:
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise PreparationError(f"refusing to overwrite output: {output}")
        if args.time_unit_ma <= 0:
            raise PreparationError("time unit must be positive")

        fasta_path = regular(args.fasta)
        topology_path = regular(args.topology)
        freeze_path = regular(args.topology_freeze)
        gate_path = regular(args.dating_gate)
        constraint_path = regular(args.constraints)
        freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        if freeze.get("status") != "PASS":
            raise PreparationError("topology freeze is not PASS")
        if freeze.get("accepted_topology", {}).get("sha256") != sha256(topology_path):
            raise PreparationError("topology does not match the frozen topology")
        if gate.get("status") != "PASS_SECONDARY_TIMETREE_CALIBRATION":
            raise PreparationError("dating gate is not PASS_SECONDARY_TIMETREE_CALIBRATION")
        if gate.get("production_mcmctree_allowed") is not True:
            raise PreparationError("dating gate does not allow MCMCTree")
        if gate.get("secondary_constraints", {}).get("sha256") != sha256(constraint_path):
            raise PreparationError("constraints do not match the dating-gate binding")

        records = read_fasta(fasta_path)
        tree = Phylo.read(str(topology_path), "newick")
        tip_order = [tip.name for tip in tree.get_terminals()]
        if any(not tip for tip in tip_order) or len(tip_order) != len(set(tip_order)):
            raise PreparationError("tree tips are missing or duplicated")
        if set(tip_order) != set(records):
            raise PreparationError("frozen topology and codon supermatrix tips differ")
        tip_by_name = {tip.name: tip for tip in tree.get_terminals()}

        rows = [row for row in read_tsv(constraint_path) if active_status(row.get("status", ""))]
        gate_ids = set(gate.get("active_secondary_constraint_ids", []))
        if not rows or {row["constraint_id"] for row in rows} != gate_ids:
            raise PreparationError("active constraints differ from the dating gate")
        calibration_by_clade: dict[int, tuple[float, float]] = {}
        calibration_manifest: list[dict[str, object]] = []
        clade_tip_sets: dict[int, set[str]] = {}
        for row in rows:
            descendant_a = row["descendant_a"]
            descendant_b = row["descendant_b"]
            if descendant_a not in tip_by_name or descendant_b not in tip_by_name:
                raise PreparationError(f"constraint tips not found: {row['constraint_id']}")
            clade = tree.common_ancestor(tip_by_name[descendant_a], tip_by_name[descendant_b])
            if id(clade) in calibration_by_clade:
                raise PreparationError(f"multiple constraints target one MRCA: {row['constraint_id']}")
            minimum_ma = float(row["minimum_ma"])
            maximum_ma = float(row["maximum_ma"])
            if not 0 < minimum_ma < maximum_ma:
                raise PreparationError(f"invalid calibration interval: {row['constraint_id']}")
            scaled = (minimum_ma / args.time_unit_ma, maximum_ma / args.time_unit_ma)
            calibration_by_clade[id(clade)] = scaled
            clade_tip_sets[id(clade)] = {tip.name for tip in clade.get_terminals()}
            calibration_manifest.append(
                {
                    "constraint_id": row["constraint_id"],
                    "node_label": row["node_label"],
                    "descendant_a": descendant_a,
                    "descendant_b": descendant_b,
                    "minimum_ma": minimum_ma,
                    "maximum_ma": maximum_ma,
                    "minimum_mcmctree_units": scaled[0],
                    "maximum_mcmctree_units": scaled[1],
                    "mrca_tip_count": len(clade_tip_sets[id(clade)]),
                }
            )
        for ancestor_id, ancestor_interval in calibration_by_clade.items():
            for descendant_id, descendant_interval in calibration_by_clade.items():
                if ancestor_id == descendant_id:
                    continue
                ancestor_tips = clade_tip_sets[ancestor_id]
                descendant_tips = clade_tip_sets[descendant_id]
                if descendant_tips < ancestor_tips and ancestor_interval[1] <= descendant_interval[0]:
                    raise PreparationError("nested calibration intervals have no chronological solution")

        version, executable_sha256 = probe_mcmctree(args.mcmctree)
        if version != "4.10.10":
            raise PreparationError(f"unexpected MCMCTree version: {version}")
        root_upper = max(maximum for minimum, maximum in calibration_by_clade.values()) * 1.20
        alignment_length = len(next(iter(records.values())))

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            phylip_path = staging / "codon_positions.phy"
            partition_lengths = write_partitioned_phylip(phylip_path, records, tip_order)
            tree_path = staging / "calibrated_topology.trees"
            tree_path.write_text(
                f" {len(tip_order)} 1\n\n{render_calibrated_tree(tree, calibration_by_clade)}\n",
                encoding="ascii",
            )
            controls = {
                "prior.ctl": control_text(
                    seed=2026072000, usedata=0, outfile="prior.out", mcmcfile="prior.mcmc.txt",
                    burnin=20000, sampfreq=10, nsample=20000, root_age_upper=root_upper
                ),
                "hessian.ctl": control_text(
                    seed=2026072003, usedata=3, outfile="hessian.out", mcmcfile="hessian.mcmc.txt",
                    burnin=2000, sampfreq=10, nsample=2000, root_age_upper=root_upper
                ),
                "posterior_chain1.ctl": control_text(
                    seed=2026072001, usedata=2, outfile="posterior_chain1.out",
                    mcmcfile="posterior_chain1.mcmc.txt", burnin=100000, sampfreq=10,
                    nsample=50000, root_age_upper=root_upper
                ),
                "posterior_chain2.ctl": control_text(
                    seed=2026072002, usedata=2, outfile="posterior_chain2.out",
                    mcmcfile="posterior_chain2.mcmc.txt", burnin=100000, sampfreq=10,
                    nsample=50000, root_age_upper=root_upper
                ),
            }
            for filename, content in controls.items():
                (staging / filename).write_text(content, encoding="ascii")

            manifest = {
                "schema_version": 1,
                "workflow": "mcmctree_timetree_secondary_bundle",
                "status": "PASS_PREPARED",
                "calibration_claim": "TimeTree secondary-calibrated; not fossil-calibrated",
                "fasta": binding(fasta_path),
                "topology": binding(topology_path),
                "topology_freeze": binding(freeze_path),
                "dating_gate": binding(gate_path),
                "constraints": binding(constraint_path),
                "tip_count": len(tip_order),
                "tip_order": tip_order,
                "alignment_length_bp": alignment_length,
                "ndata": 3,
                "partition_policy": "codon_positions_1_2_3",
                "partition_lengths_bp": partition_lengths,
                "time_unit_ma": args.time_unit_ma,
                "root_age_safe_upper_mcmctree_units": root_upper,
                "active_constraint_count": len(calibration_manifest),
                "active_constraints": calibration_manifest,
                "model": "HKY85+Gamma5",
                "clock": "independent_rates",
                "mcmctree": {
                    "version": version,
                    "sha256": executable_sha256,
                    "basename": regular(args.mcmctree).name,
                },
                "run_design": {
                    "prior_only": "usedata=0",
                    "hessian": "usedata=3 creates out.BV",
                    "posterior_chain1": "usedata=2 consumes checksum-bound in.BV",
                    "posterior_chain2": "usedata=2 consumes the same checksum-bound in.BV",
                },
                "posterior_input_filename": "in.BV",
            }
            manifest_path = staging / "run_manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            checksum_rows = []
            for path in sorted(staging.iterdir(), key=lambda item: item.name):
                if path.is_file():
                    checksum_rows.append((path.name, sha256(path)))
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n" + "".join(f"{name}\t{digest}\n" for name, digest in checksum_rows),
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS_PREPARED\t{output}")
        return 0
    except (OSError, UnicodeError, ValueError, KeyError, json.JSONDecodeError, PreparationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
