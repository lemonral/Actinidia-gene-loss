#!/usr/bin/env python3
"""Download public assets from a TSV manifest and verify their integrity.

The downloader is fail-closed: destinations must remain below ``--data-root``;
enabled rows and paths must be unique; declared sizes/checksums must match; and
a JSON report is written even when transfer or integrity validation fails.

Direct downloads use curl unless a stable host is explicitly approved for
bounded aria2 segmentation. Selected foreign hosts can instead use a user-local
HTTP proxy; proxy settings are applied per request and are never exported
machine-wide. Proxy-routed rows use curl because short-lived signed redirects
are safer to refresh through the original repository URL.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from urllib.parse import urlparse


REQUIRED_COLUMNS = {
    "asset_id",
    "assembly_unit_id",
    "asset_type",
    "url",
    "relative_path",
    "expected_bytes",
    "md5",
    "sha256",
    "download",
    "source_note",
}
PROXY_ENVIRONMENT = {
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "ALL_PROXY",
    "NO_PROXY",
    "http_proxy",
    "https_proxy",
    "all_proxy",
    "no_proxy",
}


class ManifestError(RuntimeError):
    """Raised when a download row is unsafe or internally inconsistent."""


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def clean_environment() -> dict[str, str]:
    environment = os.environ.copy()
    for variable in PROXY_ENVIRONMENT:
        environment.pop(variable, None)
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[variable] = "1"
    return environment


def is_under(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def destination_path(data_root: Path, relative_path: str) -> Path:
    raw = Path(relative_path)
    if raw.is_absolute() or not relative_path.strip():
        raise ManifestError(f"relative_path must be a non-empty relative path: {relative_path!r}")
    destination = (data_root / raw).resolve()
    if not is_under(destination, data_root):
        raise ManifestError(f"relative_path escapes --data-root: {relative_path!r}")
    return destination


def host_matches(host: str, declared: list[str]) -> bool:
    return any(host == suffix or host.endswith("." + suffix) for suffix in declared)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise ManifestError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = []
        asset_ids: set[str] = set()
        relative_paths: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            asset_id = row["asset_id"]
            if not asset_id:
                raise ManifestError(f"{path}:{line_number}: empty asset_id")
            if asset_id in asset_ids:
                raise ManifestError(f"{path}:{line_number}: duplicate asset_id {asset_id!r}")
            asset_ids.add(asset_id)
            relative = row["relative_path"]
            if relative in relative_paths:
                raise ManifestError(f"{path}:{line_number}: duplicate relative_path {relative!r}")
            relative_paths.add(relative)
            if row["download"] not in {"true", "false"}:
                raise ManifestError(f"{path}:{line_number}: download must be true or false")
            parsed = urlparse(row["url"])
            if parsed.scheme not in {"http", "https"} or not parsed.hostname:
                raise ManifestError(f"{path}:{line_number}: unsupported URL {row['url']!r}")
            for field in ("md5", "sha256"):
                value = row[field]
                expected_length = 32 if field == "md5" else 64
                if value and (len(value) != expected_length or any(c not in "0123456789abcdefABCDEF" for c in value)):
                    raise ManifestError(f"{path}:{line_number}: invalid {field}")
            if row["expected_bytes"]:
                try:
                    size = int(row["expected_bytes"])
                except ValueError as exc:
                    raise ManifestError(f"{path}:{line_number}: expected_bytes is not an integer") from exc
                if size <= 0:
                    raise ManifestError(f"{path}:{line_number}: expected_bytes must be positive")
            rows.append(row)
    return rows


def verify(path: Path, row: dict[str, str], route: str) -> dict[str, object]:
    observed_size = path.stat().st_size
    expected_size = int(row["expected_bytes"]) if row["expected_bytes"] else None
    observed_md5 = digest(path, "md5") if row["md5"] else ""
    observed_sha256 = digest(path, "sha256")
    problems: list[str] = []
    if expected_size is not None and observed_size != expected_size:
        problems.append(f"size:{observed_size}!={expected_size}")
    if row["md5"] and observed_md5.lower() != row["md5"].lower():
        problems.append("md5_mismatch")
    if row["sha256"] and observed_sha256.lower() != row["sha256"].lower():
        problems.append("sha256_mismatch")
    return {
        "asset_id": row["asset_id"],
        "assembly_unit_id": row["assembly_unit_id"],
        "asset_type": row["asset_type"],
        "relative_path": row["relative_path"],
        "bytes": observed_size,
        "md5": observed_md5,
        "sha256": observed_sha256,
        "publisher_checksum_declared": bool(row["md5"] or row["sha256"]),
        "verified": not problems and bool(row["md5"] or row["sha256"]),
        "route": route,
        "status": "verified" if not problems else "integrity_failed",
        "problems": problems,
    }


def write_report(path: Path, reports: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def curl_download(
    url: str,
    part: Path,
    *,
    proxy: str | None,
    direct_domains: list[str],
    retries: int,
    environment: dict[str, str],
) -> None:
    command = [
        "curl",
        "--location",
        "--fail",
        "--show-error",
        "--retry",
        str(retries),
        "--retry-all-errors",
        "--retry-delay",
        "3",
        "--connect-timeout",
        "30",
        "--continue-at",
        "-",
        "--output",
        str(part),
    ]
    if proxy:
        command.extend(["--proxy", proxy])
        if direct_domains:
            command.extend(["--noproxy", ",".join(direct_domains)])
    command.append(url)
    subprocess.run(command, check=True, env=environment)


def aria2_download(
    url: str,
    part: Path,
    *,
    connections: int,
    proxy: str | None,
    environment: dict[str, str],
) -> None:
    executable = shutil.which("aria2c")
    if not executable:
        raise ManifestError("--connections > 1 requires aria2c on PATH")
    command = [
        executable,
        "--allow-overwrite=false",
        "--auto-file-renaming=false",
        "--continue=true",
        "--file-allocation=none",
        "--max-connection-per-server",
        str(connections),
        "--split",
        str(connections),
        "--min-split-size=1M",
        "--max-tries=6",
        "--retry-wait=5",
        "--dir",
        str(part.parent),
        "--out",
        part.name,
    ]
    if proxy:
        command.extend(["--all-proxy", proxy])
    command.append(url)
    subprocess.run(command, check=True, env=environment)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--connections", type=int, default=1, help="Direct-download segments, 1-10")
    parser.add_argument(
        "--curl-retries",
        type=int,
        default=50,
        help="Per-file curl retries for intermittent repository disconnects, 1-100",
    )
    parser.add_argument("--proxy", help="HTTP proxy used only for --proxy-domain hosts")
    parser.add_argument("--proxy-domain", action="append", default=[])
    parser.add_argument(
        "--segmented-proxy-domain",
        action="append",
        default=[],
        help=(
            "Stable proxy-routed hosts that may use aria2 segmentation and resume. "
            "Signed repository redirects should remain curl-routed."
        ),
    )
    parser.add_argument(
        "--segmented-direct-domain",
        action="append",
        default=[],
        help="Stable direct hosts that may use bounded aria2 segmentation and resume.",
    )
    parser.add_argument("--direct-domain", action="append", default=[])
    args = parser.parse_args()

    if not 1 <= args.connections <= 10:
        raise SystemExit("ERROR: --connections must be between 1 and 10")
    if not 1 <= args.curl_retries <= 100:
        raise SystemExit("ERROR: --curl-retries must be between 1 and 100")
    if args.proxy_domain and not args.proxy:
        raise SystemExit("ERROR: --proxy-domain requires --proxy")
    if args.segmented_proxy_domain and not args.proxy:
        raise SystemExit("ERROR: --segmented-proxy-domain requires --proxy")
    undeclared_segmented = [
        host for host in args.segmented_proxy_domain if not host_matches(host, args.proxy_domain)
    ]
    if undeclared_segmented:
        raise SystemExit(
            "ERROR: every --segmented-proxy-domain must also be declared with --proxy-domain: "
            + ",".join(undeclared_segmented)
        )

    data_root = args.data_root.expanduser().resolve()
    data_root.mkdir(parents=True, exist_ok=True)
    reports: list[dict[str, object]] = []
    try:
        rows = read_rows(args.manifest.expanduser().resolve())
        environment = clean_environment()
        for row in rows:
            if row["download"] != "true":
                continue
            destination = destination_path(data_root, row["relative_path"])
            destination.parent.mkdir(parents=True, exist_ok=True)
            part = destination.with_name(destination.name + ".part")
            hostname = urlparse(row["url"]).hostname or ""
            proxied = host_matches(hostname, args.proxy_domain)
            segmented_proxy = proxied and host_matches(hostname, args.segmented_proxy_domain)
            segmented_direct = (not proxied) and host_matches(
                hostname, args.segmented_direct_domain
            )
            aria2_control = part.with_name(part.name + ".aria2")
            route = "proxy" if proxied else "direct"
            if not destination.exists():
                if segmented_proxy and args.connections > 1:
                    aria2_download(
                        row["url"],
                        part,
                        connections=args.connections,
                        proxy=args.proxy,
                        environment=environment,
                    )
                elif proxied:
                    curl_download(
                        row["url"],
                        part,
                        proxy=args.proxy,
                        direct_domains=args.direct_domain,
                        retries=args.curl_retries,
                        environment=environment,
                    )
                elif segmented_direct and args.connections > 1:
                    aria2_download(
                        row["url"],
                        part,
                        connections=args.connections,
                        proxy=None,
                        environment=environment,
                    )
                elif aria2_control.exists():
                    # A segmented partial file can contain sparse, out-of-order
                    # ranges. Resume it with its aria2 control file rather than
                    # treating the apparent file size as a contiguous curl
                    # offset. One connection is deliberately used here because
                    # this branch is also the safe fallback after a host rejects
                    # multi-connection TLS sessions.
                    aria2_download(
                        row["url"],
                        part,
                        connections=1,
                        proxy=None,
                        environment=environment,
                    )
                else:
                    curl_download(
                        row["url"],
                        part,
                        proxy=None,
                        direct_domains=[],
                        retries=args.curl_retries,
                        environment=environment,
                    )
                part.replace(destination)
            report = verify(destination, row, route)
            reports.append(report)
            write_report(args.report, reports)
            if report["problems"]:
                raise ManifestError(
                    f"integrity check failed for {row['asset_id']}: "
                    f"{','.join(str(value) for value in report['problems'])}"
                )
    except (OSError, csv.Error, ManifestError, subprocess.CalledProcessError) as exc:
        reports.append({"status": "failed", "error": str(exc)})
        write_report(args.report, reports)
        raise SystemExit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
