#!/usr/bin/env python3
"""Compare diploid and polyploid gene-loss rates using all 23 genomes.

Each haplotype or subgenome is one observation. The script reports an exact
two-sided Mann-Whitney test based on all group-label assignments and a separate
exact permutation test for the difference in group means.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from statistics import mean, median


REQUIRED_COLUMNS = {
    "assembly_unit_id",
    "species",
    "ploidy",
    "ploidy_class",
    "retained",
    "decayed",
    "deleted",
    "positive_loss",
    "resolved_denominator",
    "positive_loss_rate",
}


class PloidyComparisonError(ValueError):
    """Raised when the input cannot support the declared comparison."""


@dataclass(frozen=True)
class Observation:
    assembly_unit_id: str
    species: str
    ploidy: str
    ploidy_group: str
    retained: int
    decayed: int
    deleted: int
    positive_loss: int
    resolved_denominator: int
    positive_loss_rate: float


@dataclass(frozen=True)
class ExactTests:
    mann_whitney_u_polyploid: float
    mann_whitney_u_complement: float
    mann_whitney_u_reported: float
    mann_whitney_exact_p: float
    mean_difference: float
    permutation_exact_p: float
    permutation_assignments: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-genomes", type=int, default=23)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_observations(path: Path, expected_genomes: int) -> list[Observation]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        if not REQUIRED_COLUMNS.issubset(fields):
            missing = sorted(REQUIRED_COLUMNS - fields)
            raise PloidyComparisonError(f"missing columns: {', '.join(missing)}")
        rows = [dict(row) for row in reader]

    if len(rows) != expected_genomes:
        raise PloidyComparisonError(
            f"expected {expected_genomes} genomes, found {len(rows)}"
        )

    observations: list[Observation] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        unit = row["assembly_unit_id"]
        if not unit or unit in seen:
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: missing or duplicate genome ID"
            )
        seen.add(unit)
        try:
            retained = int(row["retained"])
            decayed = int(row["decayed"])
            deleted = int(row["deleted"])
            positive = int(row["positive_loss"])
            denominator = int(row["resolved_denominator"])
            reported_rate = float(row["positive_loss_rate"])
        except ValueError as error:
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: invalid numeric value"
            ) from error
        if positive != decayed + deleted:
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: positive loss is not decayed + deleted"
            )
        if denominator != retained + decayed + deleted or denominator <= 0:
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: invalid resolved denominator"
            )
        calculated_rate = positive / denominator
        if not math.isclose(calculated_rate, reported_rate, rel_tol=0, abs_tol=1e-12):
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: reported loss rate does not close"
            )
        ploidy_class = row["ploidy_class"].strip().lower()
        if ploidy_class == "diploid":
            group = "diploid"
        elif ploidy_class in {"tetraploid", "hexaploid", "polyploid"}:
            group = "polyploid"
        else:
            raise PloidyComparisonError(
                f"{path.name}:{line_number}: unsupported ploidy class {ploidy_class!r}"
            )
        observations.append(
            Observation(
                assembly_unit_id=unit,
                species=row["species"],
                ploidy=row["ploidy"],
                ploidy_group=group,
                retained=retained,
                decayed=decayed,
                deleted=deleted,
                positive_loss=positive,
                resolved_denominator=denominator,
                positive_loss_rate=calculated_rate,
            )
        )
    return observations


def average_ranks(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda pair: pair[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        average_rank = (index + 1 + end) / 2
        for tied in range(index, end):
            ranks[ordered[tied][0]] = average_rank
        index = end
    return ranks


def exact_tests(polyploid: list[float], diploid: list[float]) -> ExactTests:
    combined = polyploid + diploid
    n_polyploid = len(polyploid)
    n_diploid = len(diploid)
    if not polyploid or not diploid:
        raise PloidyComparisonError("both ploidy groups must be nonempty")
    ranks = average_ranks(combined)
    rank_sum = sum(ranks[:n_polyploid])
    u_polyploid = rank_sum - n_polyploid * (n_polyploid + 1) / 2
    u_complement = n_polyploid * n_diploid - u_polyploid
    observed_minimum_u = min(u_polyploid, u_complement)
    observed_mean_difference = mean(polyploid) - mean(diploid)
    combined_sum = sum(combined)

    assignments = math.comb(len(combined), n_polyploid)
    mann_whitney_extreme = 0
    mean_difference_extreme = 0
    for selected in combinations(range(len(combined)), n_polyploid):
        selected_rank_sum = sum(ranks[index] for index in selected)
        selected_u = selected_rank_sum - n_polyploid * (n_polyploid + 1) / 2
        selected_minimum_u = min(
            selected_u,
            n_polyploid * n_diploid - selected_u,
        )
        if selected_minimum_u <= observed_minimum_u + 1e-12:
            mann_whitney_extreme += 1

        selected_sum = sum(combined[index] for index in selected)
        permuted_difference = (
            selected_sum / n_polyploid
            - (combined_sum - selected_sum) / n_diploid
        )
        if abs(permuted_difference) + 1e-15 >= abs(observed_mean_difference):
            mean_difference_extreme += 1

    return ExactTests(
        mann_whitney_u_polyploid=u_polyploid,
        mann_whitney_u_complement=u_complement,
        mann_whitney_u_reported=observed_minimum_u,
        mann_whitney_exact_p=mann_whitney_extreme / assignments,
        mean_difference=observed_mean_difference,
        permutation_exact_p=mean_difference_extreme / assignments,
        permutation_assignments=assignments,
    )


def scenario_rows(observations: list[Observation]) -> list[dict[str, object]]:
    scenarios = (
        ("all_genomes", False, False),
        ("exclude_zhejiangensis", True, False),
        ("exclude_macrosperma", False, True),
        ("exclude_zhejiangensis_and_macrosperma", True, True),
    )
    output: list[dict[str, object]] = []
    for name, exclude_zhejiangensis, exclude_macrosperma in scenarios:
        selected = [
            row
            for row in observations
            if not (
                exclude_zhejiangensis
                and "zhejiangensis" in row.species.lower()
            )
            and not (
                exclude_macrosperma
                and row.species == "Actinidia macrosperma"
            )
        ]
        polyploid = [
            row.positive_loss_rate
            for row in selected
            if row.ploidy_group == "polyploid"
        ]
        diploid = [
            row.positive_loss_rate
            for row in selected
            if row.ploidy_group == "diploid"
        ]
        tests = exact_tests(polyploid, diploid)
        output.append(
            {
                "analysis": name,
                "polyploid_n": len(polyploid),
                "diploid_n": len(diploid),
                "polyploid_mean": mean(polyploid),
                "diploid_mean": mean(diploid),
                "polyploid_median": median(polyploid),
                "diploid_median": median(diploid),
                "mean_difference_percentage_points": tests.mean_difference * 100,
                **asdict(tests),
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    observations = read_observations(args.unit_summary, args.expected_genomes)
    summaries = scenario_rows(observations)
    primary = summaries[0]
    if primary["polyploid_n"] != 11 or primary["diploid_n"] != 12:
        raise PloidyComparisonError("expected 11 polyploid and 12 diploid genomes")
    if not math.isclose(primary["mann_whitney_u_reported"], 22.0, abs_tol=1e-12):
        raise PloidyComparisonError("primary Mann-Whitney U does not match 22")
    if not math.isclose(
        primary["mann_whitney_exact_p"],
        0.00562245669258726,
        rel_tol=0,
        abs_tol=1e-15,
    ):
        raise PloidyComparisonError("primary exact Mann-Whitney P value changed")

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=args.output_dir.parent,
        prefix=f".{args.output_dir.name}.",
    ) as temporary_name:
        temporary = Path(temporary_name)
        observation_rows = [asdict(row) for row in observations]
        write_tsv(temporary / "ploidy_group_observations.tsv", observation_rows)
        write_tsv(temporary / "ploidy_group_statistics.tsv", summaries)
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_EXACT_PLOIDY_GROUP_COMPARISON",
            "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "input": {
                "basename": args.unit_summary.name,
                "sha256": sha256(args.unit_summary),
            },
            "definitions": {
                "observation": "one haplotype, subgenome, or single-genome assembly",
                "positive_loss": "decayed + deleted",
                "denominator": "retained + decayed + deleted",
                "mann_whitney": "two-sided exact label-enumeration test using the minimum U tail",
                "permutation": "two-sided exact label-enumeration test for the difference in group means",
            },
            "counts": {
                "genomes": len(observations),
                "polyploid_genomes": primary["polyploid_n"],
                "diploid_genomes": primary["diploid_n"],
            },
            "primary_statistics": primary,
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if args.output_dir.exists():
            raise PloidyComparisonError(
                f"refusing to replace existing output directory: {args.output_dir}"
            )
        temporary.rename(args.output_dir)

    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
