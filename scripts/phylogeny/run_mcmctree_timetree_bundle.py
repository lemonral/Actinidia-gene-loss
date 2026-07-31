#!/usr/bin/env python3
"""Run one checksum-bound TimeTree-secondary MCMCTree workflow sequentially."""

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
from datetime import datetime, timezone
from pathlib import Path


class RunError(RuntimeError):
    pass


STAGES = (
    ("prior", "prior.ctl"),
    ("hessian", "hessian.ctl"),
    ("posterior_chain1", "posterior_chain1.ctl"),
    ("posterior_chain2", "posterior_chain2.ctl"),
)
SHARED_INPUTS = ("codon_positions.phy", "calibrated_topology.trees")


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise RunError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def read_bundle_checksums(bundle: Path) -> dict[str, str]:
    checksum_path = regular(bundle / "checksums.tsv")
    with checksum_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise RunError("invalid prepared-bundle checksum header")
        rows = list(reader)
    observed: dict[str, str] = {}
    for row in rows:
        name = row["file"]
        digest = row["sha256"]
        if not name or Path(name).name != name or name in observed:
            raise RunError(f"invalid or duplicate checksum filename: {name!r}")
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise RunError(f"invalid checksum for {name}")
        observed[name] = digest
    expected_files = {
        path.name for path in bundle.iterdir()
        if path.is_file() and path.name != "checksums.tsv"
    }
    if set(observed) != expected_files:
        raise RunError("prepared-bundle checksum inventory does not close")
    for name, expected in observed.items():
        if sha256(regular(bundle / name)) != expected:
            raise RunError(f"prepared-bundle checksum mismatch: {name}")
    return observed


def validate_bundle(bundle: Path, executable: Path) -> tuple[dict[str, object], dict[str, str]]:
    if not bundle.is_dir() or bundle.is_symlink():
        raise RunError(f"invalid bundle directory: {bundle}")
    manifest_path = regular(bundle / "run_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_PREPARED":
        raise RunError("prepared bundle is not PASS_PREPARED")
    if manifest.get("workflow") != "mcmctree_timetree_secondary_bundle":
        raise RunError("unexpected prepared-bundle workflow")
    if manifest.get("calibration_claim") != "TimeTree secondary-calibrated; not fossil-calibrated":
        raise RunError("calibration claim is missing or changed")
    if manifest.get("active_constraint_count") != 4:
        raise RunError("expected exactly four active TimeTree constraints")
    if manifest.get("mcmctree", {}).get("version") != "4.10.10":
        raise RunError("prepared bundle is not bound to MCMCTree 4.10.10")
    if manifest.get("mcmctree", {}).get("sha256") != sha256(executable):
        raise RunError("MCMCTree executable checksum differs from prepared bundle")
    checksums = read_bundle_checksums(bundle)
    required = {"run_manifest.json", *SHARED_INPUTS, *(control for _, control in STAGES)}
    if not required <= set(checksums):
        raise RunError("prepared bundle is missing required checksum-bound inputs")
    return manifest, checksums


def probe_version(executable: Path, program: str) -> str:
    completed = subprocess.run(
        [str(executable)], input="\n", text=True, capture_output=True, check=False, timeout=30
    )
    match = re.search(
        rf"{re.escape(program)} in paml version ([0-9.]+)",
        completed.stdout + "\n" + completed.stderr,
        flags=re.IGNORECASE,
    )
    if not match:
        raise RunError(f"could not recover the {program} version banner")
    return match.group(1)


def control_outputs(control: Path) -> tuple[int, str, str]:
    values: dict[str, str] = {}
    for raw in control.read_text(encoding="ascii").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        values[key.strip()] = value.strip().split()[0]
    try:
        return int(values["usedata"]), values["outfile"], values["mcmcfile"]
    except (KeyError, ValueError) as error:
        raise RunError(f"invalid MCMCTree control file: {control}") from error


def run_stage(
    *, stage: str, control_name: str, bundle: Path, output: Path, executable: Path,
    hessian_sha256: str | None, environment: dict[str, str],
) -> dict[str, object]:
    run_dir = output / "runs" / stage
    run_dir.mkdir(parents=True)
    for name in (*SHARED_INPUTS, control_name):
        shutil.copyfile(bundle / name, run_dir / name)
    if stage.startswith("posterior_"):
        if hessian_sha256 is None:
            raise RunError("posterior stage lacks the frozen Hessian checksum")
        shutil.copyfile(output / "runs" / "hessian" / "out.BV", run_dir / "in.BV")
        if sha256(run_dir / "in.BV") != hessian_sha256:
            raise RunError(f"{stage}: in.BV differs from frozen Hessian")

    control = run_dir / control_name
    usedata, outfile_name, mcmcfile_name = control_outputs(control)
    expected_mode = {"prior": 0, "hessian": 3, "posterior_chain1": 2, "posterior_chain2": 2}[stage]
    if usedata != expected_mode:
        raise RunError(f"{stage}: unexpected usedata={usedata}")
    command = [str(executable), control_name]
    stdout = run_dir / "console.stdout"
    stderr = run_dir / "console.stderr"
    started = now()
    with stdout.open("w", encoding="utf-8") as stdout_handle, stderr.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            command,
            cwd=run_dir,
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            env=environment,
            check=False,
        )
    if completed.returncode != 0:
        raise RunError(f"{stage}: MCMCTree exited {completed.returncode}")

    required_outputs = [run_dir / outfile_name]
    if stage != "hessian":
        required_outputs.append(run_dir / mcmcfile_name)
    if stage == "hessian":
        required_outputs.append(run_dir / "out.BV")
    for path in required_outputs:
        regular(path)
    return {
        "stage": stage,
        "started_at_utc": started,
        "finished_at_utc": now(),
        "returncode": completed.returncode,
        "control_sha256": sha256(bundle / control_name),
        "stdout_sha256": sha256(stdout),
        "stderr_sha256": sha256(stderr),
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in required_outputs
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--mcmctree", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--resume-failed-missing-baseml",
        action="store_true",
        help="Resume only the registered zero-byte out.BV failure caused by BASEML missing from PATH",
    )
    parser.add_argument(
        "--reuse-prior-hessian-from",
        type=Path,
        help="Reuse exact completed prior/Hessian stages from a stopped older prepared bundle",
    )
    args = parser.parse_args()

    output = args.output_dir.expanduser().resolve()
    state_path = output / "state.json"
    state: dict[str, object] | None = None
    try:
        bundle = args.bundle.expanduser().resolve()
        executable = regular(args.mcmctree)
        if not os.access(executable, os.X_OK):
            raise RunError(f"MCMCTree is not executable: {executable}")
        baseml = regular(executable.parent / "baseml")
        if not os.access(baseml, os.X_OK):
            raise RunError(f"BASEML is not executable: {baseml}")
        version = probe_version(executable, "MCMCTREE")
        if version != "4.10.10":
            raise RunError(f"unexpected MCMCTree version: {version}")
        baseml_version = probe_version(baseml, "BASEML")
        if baseml_version != "4.10.10":
            raise RunError(f"unexpected BASEML version: {baseml_version}")
        manifest, bundle_checksums = validate_bundle(bundle, executable)
        software = {
            "mcmctree": {
                "path": str(executable), "version": version, "sha256": sha256(executable)
            },
            "baseml": {
                "path": str(baseml), "version": baseml_version, "sha256": sha256(baseml)
            },
        }
        prepared_binding = {
            "path": str(bundle),
            "manifest_sha256": bundle_checksums["run_manifest.json"],
            "checksums_sha256": sha256(bundle / "checksums.tsv"),
        }
        environment = os.environ.copy()
        environment["PATH"] = str(executable.parent) + os.pathsep + environment.get("PATH", "")

        start_index = 0
        hessian_sha256: str | None = None
        completed_stages: list[dict[str, object]] = []
        if output.exists():
            if args.reuse_prior_hessian_from is not None:
                raise RunError("cannot reuse stages into an existing output directory")
            if not args.resume_failed_missing_baseml:
                raise RunError(f"refusing to overwrite output directory: {output}")
            state = json.loads(regular(state_path).read_text(encoding="utf-8"))
            if state.get("workflow") != "sequential_mcmctree_timetree_secondary":
                raise RunError("cannot resume an unrelated workflow")
            if state.get("status") != "ERROR" or state.get("active_stage") != "hessian":
                raise RunError("resume is allowed only for the registered failed Hessian stage")
            if "missing, empty, or symlink file" not in str(state.get("error", "")) or not str(
                state.get("error", "")
            ).endswith("/runs/hessian/out.BV"):
                raise RunError("resume error is not the registered missing out.BV failure")
            if state.get("prepared_bundle") != prepared_binding:
                raise RunError("resume prepared-bundle binding differs")
            if state.get("mcmctree") != software["mcmctree"]:
                raise RunError("resume MCMCTree binding differs")
            completed_stages = list(state.get("completed", []))
            if len(completed_stages) != 1 or completed_stages[0].get("stage") != "prior":
                raise RunError("resume would repeat or skip a completed MCMCTree stage")
            for name, bound in completed_stages[0].get("outputs", {}).items():
                path = regular(output / "runs" / "prior" / name)
                if sha256(path) != bound.get("sha256") or path.stat().st_size != bound.get("bytes"):
                    raise RunError(f"resume prior output binding differs: {name}")
            failed_hessian = output / "runs" / "hessian"
            failed_bv = failed_hessian / "out.BV"
            if not failed_bv.is_file() or failed_bv.stat().st_size != 0:
                raise RunError("registered failed out.BV is not an exact zero-byte file")
            if "baseml: not found" not in (failed_hessian / "console.stderr").read_text(
                encoding="utf-8"
            ):
                raise RunError("registered Hessian stderr lacks the exact BASEML PATH failure")
            archived_hessian = output / "runs" / "hessian.initial_missing_baseml"
            if archived_hessian.exists():
                raise RunError("registered failed Hessian archive already exists")
            os.rename(failed_hessian, archived_hessian)
            state.setdefault("recoveries", []).append(
                {
                    "at_utc": now(),
                    "reason": "registered_BASEML_PATH_failure_before_out_BV",
                    "archived_run": archived_hessian.name,
                    "prior_reused_without_rerun": True,
                }
            )
            state["status"] = "running"
            state["baseml"] = software["baseml"]
            state.pop("error", None)
            state.pop("finished_at_utc", None)
            start_index = 1
            write_json(state_path, state)
        else:
            if args.resume_failed_missing_baseml:
                raise RunError("resume requested but output directory does not exist")
            output.parent.mkdir(parents=True, exist_ok=True)
            output.mkdir()
            (output / "runs").mkdir()
            state = {
                "schema_version": 1,
                "workflow": "sequential_mcmctree_timetree_secondary",
                "status": "running",
                "started_at_utc": now(),
                "calibration_claim": manifest["calibration_claim"],
                "prepared_bundle": prepared_binding,
                "mcmctree": software["mcmctree"],
                "baseml": software["baseml"],
                "completed": [],
            }
            if args.reuse_prior_hessian_from is not None:
                source = args.reuse_prior_hessian_from.expanduser().resolve()
                if not source.is_dir() or source.is_symlink():
                    raise RunError(f"invalid prior/Hessian source: {source}")
                source_state_path = regular(source / "state.json")
                source_state = json.loads(source_state_path.read_text(encoding="utf-8"))
                if source_state.get("workflow") != "sequential_mcmctree_timetree_secondary":
                    raise RunError("prior/Hessian source is an unrelated workflow")
                if source_state.get("status") != "ERROR" or source_state.get(
                    "active_stage"
                ) != "posterior_chain1":
                    raise RunError("prior/Hessian source is not stopped at the first posterior chain")
                if source_state.get("error") != "posterior_chain1: MCMCTree exited 1":
                    raise RunError("prior/Hessian source has an unexpected stopping error")
                if source_state.get("mcmctree") != software["mcmctree"] or source_state.get(
                    "baseml"
                ) != software["baseml"]:
                    raise RunError("prior/Hessian source software binding differs")
                source_completed = list(source_state.get("completed", []))
                if [row.get("stage") for row in source_completed] != ["prior", "hessian"]:
                    raise RunError("prior/Hessian source completed-stage closure failed")
                source_bundle = Path(source_state["prepared_bundle"]["path"]).resolve()
                source_checksums = read_bundle_checksums(source_bundle)
                if source_state["prepared_bundle"] != {
                    "path": str(source_bundle),
                    "manifest_sha256": source_checksums["run_manifest.json"],
                    "checksums_sha256": sha256(source_bundle / "checksums.tsv"),
                }:
                    raise RunError("prior/Hessian source prepared-bundle binding differs")
                for name in (*SHARED_INPUTS, "prior.ctl", "hessian.ctl"):
                    if bundle_checksums[name] != source_checksums.get(name):
                        raise RunError(f"new bundle changes a reused-stage input: {name}")
                for row in source_completed:
                    stage = str(row["stage"])
                    for name, bound in row.get("outputs", {}).items():
                        path = regular(source / "runs" / stage / name)
                        if sha256(path) != bound.get("sha256") or path.stat().st_size != bound.get(
                            "bytes"
                        ):
                            raise RunError(f"reused stage output binding differs: {stage}/{name}")
                failed_chain = source / "runs" / "posterior_chain1"
                if "error: file name empty." not in (failed_chain / "console.stderr").read_text(
                    encoding="utf-8"
                ):
                    raise RunError("source chain lacks the exact missing in.BV filename error")
                if (failed_chain / "posterior_chain1.mcmc.txt").exists():
                    raise RunError("source chain produced samples before its registered failure")
                for stage in ("prior", "hessian"):
                    shutil.copytree(source / "runs" / stage, output / "runs" / stage)
                completed_stages = source_completed
                state["completed"] = completed_stages
                hessian_sha256 = sha256(output / "runs" / "hessian" / "out.BV")
                if hessian_sha256 != source_state.get("hessian_sha256"):
                    raise RunError("reused Hessian checksum differs from source state")
                state["hessian_sha256"] = hessian_sha256
                state["reused_stages"] = {
                    "source": str(source),
                    "source_state_sha256": sha256(source_state_path),
                    "stages": ["prior", "hessian"],
                    "reason": "corrected_explicit_inBV_filename_before_any_posterior_samples",
                }
                start_index = 2
            write_json(state_path, state)

        for stage, control_name in STAGES[start_index:]:
            state["active_stage"] = stage
            write_json(state_path, state)
            row = run_stage(
                stage=stage,
                control_name=control_name,
                bundle=bundle,
                output=output,
                executable=executable,
                hessian_sha256=hessian_sha256,
                environment=environment,
            )
            if stage == "hessian":
                hessian_sha256 = sha256(output / "runs" / "hessian" / "out.BV")
                state["hessian_sha256"] = hessian_sha256
            completed_stages.append(row)
            state["completed"] = completed_stages
            if read_bundle_checksums(bundle) != bundle_checksums:
                raise RunError("prepared bundle changed during the run")
            if sha256(executable) != state["mcmctree"]["sha256"]:
                raise RunError("MCMCTree executable changed during the run")
            if sha256(baseml) != state["baseml"]["sha256"]:
                raise RunError("BASEML executable changed during the run")
            write_json(state_path, state)

        state.pop("active_stage", None)
        state["status"] = "PASS_RUN_COMPLETE"
        state["finished_at_utc"] = now()
        write_json(state_path, state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError, subprocess.SubprocessError, RunError) as error:
        if state is not None:
            state["status"] = "ERROR"
            state["error"] = str(error)
            state["finished_at_utc"] = now()
            write_json(state_path, state)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
