#!/usr/bin/env python3
"""Fit species-level PGLS to lineage-specific/non-shared gene-loss rates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.io_utils import SchemaError  # noqa: E402
from geneloss_repro.pgls import parse_named_sensitivities, run_species_pgls  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path, help="One validated row per biological species")
    parser.add_argument("--time-tree", required=True, type=Path, help="Rooted ultrametric branch-length Newick")
    parser.add_argument(
        "--input-pass-report",
        required=True,
        type=Path,
        help="Checksum-bound species-PGLS input-builder PASS report",
    )
    parser.add_argument(
        "--species-loss-manifest",
        required=True,
        type=Path,
        help="Schema-2.0 PASS manifest from biological-species loss aggregation",
    )
    parser.add_argument(
        "--ploidy-ledger-pass-report",
        required=True,
        type=Path,
        help="Checksum-bound biological-species ploidy-ledger PASS report",
    )
    parser.add_argument(
        "--time-tree-pass-report",
        required=True,
        type=Path,
        help="Checksum-bound rooted biological-species time-tree PASS report",
    )
    parser.add_argument("--predictor-column", required=True, help="Single numeric predictor already on the intended scale")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--sensitivity",
        action="append",
        default=[],
        metavar="NAME=SPECIES[,SPECIES]",
        help="Repeat for named exclusion sensitivity models",
    )
    parser.add_argument("--ultrametric-tolerance", type=float, default=1e-6)
    parser.add_argument("--species-column", default="biological_species")
    parser.add_argument(
        "--count-column", default="lineage_specific_nonshared_positive_loss_count"
    )
    parser.add_argument("--denominator-column", default="callable_denominator")
    parser.add_argument("--scope-column", default="loss_scope")
    parser.add_argument("--level-column", default="analysis_level")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        sensitivities = parse_named_sensitivities(args.sensitivity)
        outputs = run_species_pgls(
            data_path=args.data,
            tree_path=args.time_tree,
            input_pass_report_path=args.input_pass_report,
            species_loss_manifest_path=args.species_loss_manifest,
            ploidy_ledger_pass_report_path=args.ploidy_ledger_pass_report,
            tree_pass_report_path=args.time_tree_pass_report,
            output_dir=args.output_dir,
            predictor_column=args.predictor_column,
            sensitivities=sensitivities,
            ultrametric_tolerance=args.ultrametric_tolerance,
            species_column=args.species_column,
            count_column=args.count_column,
            denominator_column=args.denominator_column,
            scope_column=args.scope_column,
            level_column=args.level_column,
        )
    except (OSError, SchemaError, RuntimeError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")
    for label, path in sorted(outputs.items()):
        print(f"{label}\t{path}")


if __name__ == "__main__":
    main()
