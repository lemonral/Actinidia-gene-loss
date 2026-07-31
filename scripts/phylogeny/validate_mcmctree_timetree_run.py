#!/usr/bin/env python3
"""Validate two MCMCTree chains and publish one pooled ultrametric tree in Ma."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio import Phylo


class ValidationError(RuntimeError):
    pass


CALIBRATION_CLAIM = "TimeTree secondary-calibrated; not fossil-calibrated"
STAGES = (
    ("prior", "prior.ctl"),
    ("hessian", "hessian.ctl"),
    ("posterior_chain1", "posterior_chain1.ctl"),
    ("posterior_chain2", "posterior_chain2.ctl"),
)
SHARED_INPUTS = ("codon_positions.phy", "calibrated_topology.trees")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_path(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise ValidationError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise ValidationError(f"empty file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = file_path(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_bundle_checksums(bundle: Path) -> dict[str, str]:
    checksum_path = file_path(bundle / "checksums.tsv")
    with checksum_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise ValidationError("invalid prepared-bundle checksum header")
        rows = list(reader)
    checksums: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if not name or Path(name).name != name or name in checksums:
            raise ValidationError(f"invalid or duplicate bundle filename: {name!r}")
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValidationError(f"invalid bundle checksum: {name}")
        checksums[name] = digest
    inventory = {
        path.name for path in bundle.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise ValidationError("prepared-bundle checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(file_path(bundle / name)) != digest:
            raise ValidationError(f"prepared-bundle checksum mismatch: {name}")
    return checksums


def read_control(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in file_path(path).read_text(encoding="ascii").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        values[key.strip()] = value.strip().split()[0]
    return values


def validate_run(run_dir: Path, bundle: Path) -> tuple[dict[str, object], dict[str, object]]:
    if not run_dir.is_dir() or run_dir.is_symlink():
        raise ValidationError(f"invalid run directory: {run_dir}")
    if not bundle.is_dir() or bundle.is_symlink():
        raise ValidationError(f"invalid prepared bundle: {bundle}")
    state_path = file_path(run_dir / "state.json")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    manifest_path = file_path(bundle / "run_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checksums = read_bundle_checksums(bundle)
    if state.get("status") != "PASS_RUN_COMPLETE":
        raise ValidationError("MCMCTree run is not PASS_RUN_COMPLETE")
    if state.get("workflow") != "sequential_mcmctree_timetree_secondary":
        raise ValidationError("unexpected MCMCTree run workflow")
    if state.get("calibration_claim") != CALIBRATION_CLAIM:
        raise ValidationError("run calibration claim changed")
    if manifest.get("status") != "PASS_PREPARED" or manifest.get("workflow") != (
        "mcmctree_timetree_secondary_bundle"
    ):
        raise ValidationError("prepared MCMCTree bundle is not PASS_PREPARED")
    if manifest.get("calibration_claim") != CALIBRATION_CLAIM:
        raise ValidationError("prepared-bundle calibration claim changed")
    prepared = {
        "path": str(bundle),
        "manifest_sha256": checksums.get("run_manifest.json"),
        "checksums_sha256": sha256(bundle / "checksums.tsv"),
    }
    if state.get("prepared_bundle") != prepared:
        raise ValidationError("run/prepared-bundle binding differs")
    if manifest.get("tip_count") != len(manifest.get("tip_order", [])):
        raise ValidationError("prepared-bundle tip count/order differs")
    if manifest.get("time_unit_ma") is None or float(manifest["time_unit_ma"]) <= 0:
        raise ValidationError("invalid MCMCTree time unit")

    for program in ("mcmctree", "baseml"):
        record = state.get(program, {})
        executable = file_path(Path(str(record.get("path", ""))))
        if sha256(executable) != record.get("sha256"):
            raise ValidationError(f"{program} executable checksum changed")
        if record.get("version") != "4.10.10":
            raise ValidationError(f"unexpected {program} version")

    completed = list(state.get("completed", []))
    if [row.get("stage") for row in completed] != [stage for stage, _ in STAGES]:
        raise ValidationError("completed-stage closure/order failed")
    completed_by_stage = {str(row["stage"]): row for row in completed}
    for stage, control_name in STAGES:
        stage_dir = run_dir / "runs" / stage
        if not stage_dir.is_dir() or stage_dir.is_symlink():
            raise ValidationError(f"invalid stage directory: {stage}")
        for name in (*SHARED_INPUTS, control_name):
            if sha256(file_path(stage_dir / name)) != checksums.get(name):
                raise ValidationError(f"stage input differs from bundle: {stage}/{name}")
        row = completed_by_stage[stage]
        if row.get("returncode") != 0 or row.get("control_sha256") != checksums[control_name]:
            raise ValidationError(f"invalid completed-stage record: {stage}")
        for stream in ("stdout", "stderr"):
            stream_path = file_path(stage_dir / f"console.{stream}", allow_empty=True)
            if sha256(stream_path) != row.get(f"{stream}_sha256"):
                raise ValidationError(f"stage {stream} checksum mismatch: {stage}")
        for name, expected in row.get("outputs", {}).items():
            output_path = file_path(stage_dir / name)
            if output_path.stat().st_size != expected.get("bytes") or sha256(output_path) != (
                expected.get("sha256")
            ):
                raise ValidationError(f"stage output binding mismatch: {stage}/{name}")
    hessian = file_path(run_dir / "runs" / "hessian" / "out.BV")
    if sha256(hessian) != state.get("hessian_sha256"):
        raise ValidationError("run Hessian checksum differs")
    for stage in ("posterior_chain1", "posterior_chain2"):
        if sha256(file_path(run_dir / "runs" / stage / "in.BV")) != sha256(hessian):
            raise ValidationError(f"{stage} did not consume the frozen Hessian")
    reused = state.get("reused_stages")
    if reused:
        source = Path(str(reused.get("source", ""))).resolve()
        source_state = file_path(source / "state.json")
        if sha256(source_state) != reused.get("source_state_sha256"):
            raise ValidationError("reused prior/Hessian source-state binding changed")
        if reused.get("stages") != ["prior", "hessian"]:
            raise ValidationError("unexpected reused-stage declaration")
    return state, manifest


def load_chain(path: Path, control: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    source = file_path(path)
    with source.open(encoding="ascii") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
    if len(header) < 3 or header[0] != "Gen" or len(header) != len(set(header)):
        raise ValidationError(f"invalid MCMC header: {source}")
    try:
        data = np.loadtxt(source, delimiter="\t", skiprows=1, dtype=np.float64, ndmin=2)
    except ValueError as error:
        raise ValidationError(f"invalid numeric MCMC table: {source}") from error
    if data.shape[1] != len(header) or not np.isfinite(data).all():
        raise ValidationError(f"MCMC shape/finiteness failure: {source}")
    values = read_control(control)
    try:
        nsample = int(values["nsample"])
        sampfreq = int(values["sampfreq"])
    except (KeyError, ValueError) as error:
        raise ValidationError(f"invalid posterior control: {control}") from error
    if nsample <= 0 or sampfreq <= 0 or data.shape[0] != nsample + 1:
        raise ValidationError(f"MCMC sample count differs from control: {source}")
    generations = data[:, 0]
    expected = np.concatenate(([1.0], np.arange(1, nsample + 1, dtype=np.float64) * sampfreq))
    if not np.array_equal(generations, expected):
        raise ValidationError(f"MCMC generation sequence differs from control: {source}")
    return header[1:], data[:, 1:], generations


def autocorrelation_ess(values: np.ndarray) -> float:
    series = np.asarray(values, dtype=np.float64)
    count = series.size
    centered = series - series.mean()
    variance = float(np.dot(centered, centered) / count)
    if variance == 0:
        return float(count)
    fft_length = 1 << (2 * count - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=fft_length)
    autocovariance = np.fft.irfft(spectrum * np.conjugate(spectrum), n=fft_length)[:count]
    autocovariance /= np.arange(count, 0, -1, dtype=np.float64)
    rho = autocovariance / autocovariance[0]
    paired_sum = 0.0
    previous = math.inf
    for index in range(1, count - 1, 2):
        pair = float(rho[index] + rho[index + 1])
        if not math.isfinite(pair) or pair <= 0:
            break
        pair = min(pair, previous)
        previous = pair
        paired_sum += pair
    tau = max(1.0, 1.0 + 2.0 * paired_sum)
    return min(float(count), float(count) / tau)


def split_rhat(chains: list[np.ndarray]) -> np.ndarray:
    if len(chains) < 2 or any(chain.ndim != 2 for chain in chains):
        raise ValidationError("split-Rhat requires at least two two-dimensional chains")
    half = min(chain.shape[0] // 2 for chain in chains)
    if half < 20:
        raise ValidationError("too few samples for split-Rhat")
    split = np.stack([part for chain in chains for part in (chain[:half], chain[-half:])])
    within = np.mean(np.var(split, axis=1, ddof=1), axis=0)
    between = half * np.var(np.mean(split, axis=1), axis=0, ddof=1)
    estimated = ((half - 1.0) / half) * within + between / half
    result = np.empty_like(within)
    variable = within > 0
    result[variable] = np.sqrt(estimated[variable] / within[variable])
    result[~variable] = np.where(between[~variable] == 0, 1.0, np.inf)
    return result


def parse_numbered_tree(report: Path) -> object:
    lines = file_path(report).read_text(encoding="utf-8", errors="strict").splitlines()
    marker = next(
        (index for index, line in enumerate(lines) if line.startswith("Species tree for FigTree.")),
        None,
    )
    if marker is None:
        raise ValidationError(f"missing FigTree section: {report}")
    numbered = next((line.strip() for line in lines[marker + 1 :] if line.strip().endswith(";")), None)
    if numbered is None:
        raise ValidationError(f"missing numbered topology: {report}")
    tree = Phylo.read(io.StringIO(numbered), "newick")
    for tip in tree.get_terminals():
        if not tip.name or not re.fullmatch(r"\d+_.+", tip.name):
            raise ValidationError(f"invalid numbered terminal: {tip.name!r}")
        tip.name = tip.name.split("_", 1)[1]
    for clade in tree.get_nonterminals():
        candidate: str | None = str(clade.name) if clade.name is not None else None
        if candidate is None and clade.confidence is not None:
            confidence = float(clade.confidence)
            if confidence.is_integer():
                candidate = str(int(confidence))
        if candidate is None or not candidate.isdigit():
            raise ValidationError(
                f"invalid internal node id: name={clade.name!r}, confidence={clade.confidence!r}"
            )
        clade.name = candidate
        clade.confidence = None
    return tree


def read_reported_means(report: Path) -> dict[str, float]:
    observed: dict[str, float] = {}
    pattern = re.compile(r"^(t_n\d+)\s+([0-9.eE+-]+)\s+\(")
    for line in file_path(report).read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            observed[match.group(1)] = float(match.group(2))
    if not observed:
        raise ValidationError(f"missing reported posterior node means: {report}")
    return observed


def rooted_clades(tree: object) -> set[frozenset[str]]:
    return {
        frozenset(tip.name for tip in clade.get_terminals())
        for clade in tree.get_nonterminals()
        if clade is not tree.root
    }


def parent_by_clade(tree: object) -> dict[int, object]:
    return {id(child): parent for parent in tree.find_clades() for child in parent.clades}


def render_dated_tree(tree: object) -> str:
    def render(clade: object) -> str:
        if clade.clades:
            body = "(" + ",".join(render(child) for child in clade.clades) + ")"
        else:
            if not clade.name or re.search(r"[\s(),:;]", clade.name):
                raise ValidationError(f"unsafe dated-tree tip name: {clade.name!r}")
            body = clade.name
        if clade is not tree.root:
            if clade.branch_length is None or clade.branch_length < 0:
                raise ValidationError("dated tree has a missing or negative branch")
            body += f":{clade.branch_length:.10f}"
        return body

    return render(tree.root) + ";\n"


def tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-effective-sample-size", type=float, default=200.0)
    parser.add_argument("--maximum-split-rhat", type=float, default=1.01)
    parser.add_argument("--ultrametric-tolerance-ma", type=float, default=1e-6)
    args = parser.parse_args()
    try:
        if args.minimum_effective_sample_size <= 0:
            raise ValidationError("minimum ESS must be positive")
        if args.maximum_split_rhat < 1 or args.ultrametric_tolerance_ma <= 0:
            raise ValidationError("invalid split-Rhat or ultrametric threshold")
        run_dir = args.run_dir.expanduser().resolve()
        bundle = args.bundle.expanduser().resolve()
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise ValidationError(f"refusing to overwrite output: {output}")
        state, manifest = validate_run(run_dir, bundle)

        chains: list[np.ndarray] = []
        chain_parameters: list[str] | None = None
        numbered_trees: list[object] = []
        reports: list[Path] = []
        for stage in ("posterior_chain1", "posterior_chain2"):
            stage_dir = run_dir / "runs" / stage
            parameters, samples, _ = load_chain(
                stage_dir / f"{stage}.mcmc.txt", bundle / f"{stage}.ctl"
            )
            if chain_parameters is not None and parameters != chain_parameters:
                raise ValidationError("posterior chain headers differ")
            chain_parameters = parameters
            chains.append(samples)
            report = stage_dir / f"{stage}.out"
            reports.append(report)
            numbered_trees.append(parse_numbered_tree(report))
        assert chain_parameters is not None

        tip_count = int(manifest["tip_count"])
        age_parameters = [f"t_n{node}" for node in range(tip_count + 1, 2 * tip_count)]
        if chain_parameters[: len(age_parameters)] != age_parameters:
            raise ValidationError("MCMC internal-node columns do not match the tip count")
        if set(manifest["tip_order"]) != {tip.name for tip in numbered_trees[0].get_terminals()}:
            raise ValidationError("MCMC numbered-tree tips differ from the prepared bundle")
        if rooted_clades(numbered_trees[0]) != rooted_clades(numbered_trees[1]):
            raise ValidationError("posterior-chain numbered topologies differ")
        node_ids = {str(clade.name) for clade in numbered_trees[0].get_nonterminals()}
        if node_ids != {str(node) for node in range(tip_count + 1, 2 * tip_count)}:
            raise ValidationError("numbered topology does not contain the exact internal-node ids")

        for report, samples in zip(reports, chains, strict=True):
            reported = read_reported_means(report)
            for index, parameter in enumerate(age_parameters):
                if parameter not in reported or abs(samples[:, index].mean() - reported[parameter]) > 6e-4:
                    raise ValidationError(f"reported/sample posterior mean differs: {parameter}")

        rhats = split_rhat(chains)
        convergence_rows: list[dict[str, object]] = []
        for index, parameter in enumerate(chain_parameters):
            chain_ess = [autocorrelation_ess(chain[:, index]) for chain in chains]
            convergence_rows.append(
                {
                    "parameter": parameter,
                    "parameter_class": "node_age" if parameter in age_parameters else "model",
                    "chain1_mean": f"{chains[0][:, index].mean():.10g}",
                    "chain2_mean": f"{chains[1][:, index].mean():.10g}",
                    "chain1_ess": f"{chain_ess[0]:.6f}",
                    "chain2_ess": f"{chain_ess[1]:.6f}",
                    "combined_ess": f"{sum(chain_ess):.6f}",
                    "split_rhat": f"{rhats[index]:.8f}",
                }
            )
        minimum_ess = min(float(row["combined_ess"]) for row in convergence_rows)
        maximum_rhat = max(float(row["split_rhat"]) for row in convergence_rows)
        if minimum_ess < args.minimum_effective_sample_size:
            raise ValidationError(
                f"minimum combined ESS {minimum_ess:.3f} is below {args.minimum_effective_sample_size}"
            )
        if maximum_rhat > args.maximum_split_rhat:
            raise ValidationError(
                f"maximum split-Rhat {maximum_rhat:.6f} exceeds {args.maximum_split_rhat}"
            )

        tree = numbered_trees[0]
        node_by_id = {str(clade.name): clade for clade in tree.get_nonterminals()}
        parents = parent_by_clade(tree)
        parameter_index = {name: index for index, name in enumerate(chain_parameters)}
        pooled = np.concatenate(chains, axis=0)
        node_age_rows: list[dict[str, object]] = []
        mean_age_units: dict[str, float] = {}
        for node_id in sorted(node_by_id, key=int):
            parameter = f"t_n{node_id}"
            index = parameter_index[parameter]
            values = pooled[:, index]
            mean_age_units[node_id] = float(values.mean())
            diagnostic = convergence_rows[index]
            node_age_rows.append(
                {
                    "node_id": f"n{node_id}",
                    "descendant_tip_count": len(node_by_id[node_id].get_terminals()),
                    "descendant_tips": ",".join(sorted(tip.name for tip in node_by_id[node_id].get_terminals())),
                    "mean_ma": f"{values.mean() * float(manifest['time_unit_ma']):.8f}",
                    "q025_ma": f"{np.quantile(values, 0.025) * float(manifest['time_unit_ma']):.8f}",
                    "q975_ma": f"{np.quantile(values, 0.975) * float(manifest['time_unit_ma']):.8f}",
                    "chain1_mean_ma": f"{chains[0][:, index].mean() * float(manifest['time_unit_ma']):.8f}",
                    "chain2_mean_ma": f"{chains[1][:, index].mean() * float(manifest['time_unit_ma']):.8f}",
                    "combined_ess": diagnostic["combined_ess"],
                    "split_rhat": diagnostic["split_rhat"],
                }
            )

        scale = float(manifest["time_unit_ma"])
        for clade in tree.find_clades(order="preorder"):
            if clade is tree.root:
                clade.branch_length = None
                continue
            parent = parents[id(clade)]
            parent_age = mean_age_units[str(parent.name)]
            child_age = mean_age_units[str(clade.name)] if clade.clades else 0.0
            branch = (parent_age - child_age) * scale
            if branch <= 0:
                raise ValidationError("pooled node means do not produce positive branches")
            clade.branch_length = branch
        distances = [tree.distance(tree.root, tip) for tip in tree.get_terminals()]
        ultrametric_deviation = max(distances) - min(distances)
        if ultrametric_deviation > args.ultrametric_tolerance_ma:
            raise ValidationError("pooled dated tree is not ultrametric within tolerance")

        calibration_rows: list[dict[str, object]] = []
        for constraint in manifest.get("active_constraints", []):
            mrca = tree.common_ancestor(constraint["descendant_a"], constraint["descendant_b"])
            node_id = str(mrca.name)
            values = pooled[:, parameter_index[f"t_n{node_id}"]] * scale
            mean = float(values.mean())
            minimum = float(constraint["minimum_ma"])
            maximum = float(constraint["maximum_ma"])
            calibration_rows.append(
                {
                    "constraint_id": constraint["constraint_id"],
                    "node_label": constraint["node_label"],
                    "node_id": f"n{node_id}",
                    "minimum_ma": f"{minimum:.10g}",
                    "maximum_ma": f"{maximum:.10g}",
                    "posterior_mean_ma": f"{mean:.8f}",
                    "posterior_q025_ma": f"{np.quantile(values, 0.025):.8f}",
                    "posterior_q975_ma": f"{np.quantile(values, 0.975):.8f}",
                    "mean_inside_secondary_interval": str(minimum <= mean <= maximum).lower(),
                }
            )

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            (staging / "dated_tree.mean_ma.tre").write_text(render_dated_tree(tree), encoding="ascii")
            convergence_fields = [
                "parameter", "parameter_class", "chain1_mean", "chain2_mean", "chain1_ess",
                "chain2_ess", "combined_ess", "split_rhat",
            ]
            tsv(staging / "convergence.tsv", convergence_fields, convergence_rows)
            node_fields = [
                "node_id", "descendant_tip_count", "descendant_tips", "mean_ma", "q025_ma",
                "q975_ma", "chain1_mean_ma", "chain2_mean_ma", "combined_ess", "split_rhat",
            ]
            tsv(staging / "node_ages_ma.tsv", node_fields, node_age_rows)
            calibration_fields = [
                "constraint_id", "node_label", "node_id", "minimum_ma", "maximum_ma",
                "posterior_mean_ma", "posterior_q025_ma", "posterior_q975_ma",
                "mean_inside_secondary_interval",
            ]
            tsv(staging / "secondary_calibration_summary.tsv", calibration_fields, calibration_rows)
            validation = {
                "schema_version": 1,
                "workflow": "mcmctree_secondary_two_chain_validation_and_ultrametric_publication",
                "status": "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC",
                "calibration_claim": CALIBRATION_CLAIM,
                "run_state": binding(run_dir / "state.json"),
                "prepared_bundle_manifest": binding(bundle / "run_manifest.json"),
                "prepared_bundle_checksums": binding(bundle / "checksums.tsv"),
                "chain_count": 2,
                "samples_per_chain_including_initial_post_burnin_state": chains[0].shape[0],
                "parameter_count": len(chain_parameters),
                "node_age_parameter_count": len(age_parameters),
                "diagnostics": {
                    "ess_method": "per-chain FFT autocorrelation with Geyer initial-positive monotone pairs; summed across chains",
                    "rhat_method": "classical four-way split-Rhat from two chains",
                    "minimum_combined_ess": minimum_ess,
                    "required_minimum_combined_ess": args.minimum_effective_sample_size,
                    "maximum_split_rhat": maximum_rhat,
                    "allowed_maximum_split_rhat": args.maximum_split_rhat,
                },
                "dated_tree": {
                    "time_unit": "Ma",
                    "construction": "pooled posterior mean internal-node ages from both chains",
                    "tip_count": tip_count,
                    "root_age_ma": mean_age_units[str(tree.root.name)] * scale,
                    "maximum_root_to_tip_deviation_ma": ultrametric_deviation,
                    "allowed_root_to_tip_deviation_ma": args.ultrametric_tolerance_ma,
                    "topology_matches_between_chains": True,
                    "tip_closure": True,
                    "positive_branch_lengths": True,
                },
                "secondary_calibrations": {
                    "count": len(calibration_rows),
                    "all_posterior_means_inside_declared_intervals": all(
                        row["mean_inside_secondary_interval"] == "true" for row in calibration_rows
                    ),
                    "interpretation": "secondary TimeTree soft bounds; not fossil calibrations",
                },
                "reused_prior_hessian_exactly_bound": bool(state.get("reused_stages")),
                "outputs": [
                    "dated_tree.mean_ma.tre", "convergence.tsv", "node_ages_ma.tsv",
                    "secondary_calibration_summary.tsv",
                ],
            }
            (staging / "validation.json").write_text(
                json.dumps(validation, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            rows = [
                (path.name, sha256(path)) for path in sorted(staging.iterdir()) if path.is_file()
            ]
            (staging / "checksums.tsv").write_text(
                "file\tsha256\n" + "".join(f"{name}\t{digest}\n" for name, digest in rows),
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"PASS_MCMCTREE_VALIDATED_ULTRAMETRIC\t{output}")
        return 0
    except (
        OSError, UnicodeError, ValueError, KeyError, json.JSONDecodeError, ValidationError
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
