"""Portable plotting helpers replacing hard-coded figure scripts.

Figures have no fixed font file, path, species list, or output directory.
The exact data table used is always copied beside the figure by the workflow
caller or recorded in the run manifest.
"""

from __future__ import annotations

from pathlib import Path

from .io_utils import SchemaError, natural_key, read_tsv


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting needs matplotlib; install with: python -m pip install -e '.[plots]'") from exc
    return plt


def plot_loss_summary(summary_path: str | Path, output_path: str | Path, title: str = "") -> Path:
    """Create a stacked count plot from ``summarize-loss`` output."""
    plt = _plt()
    rows = read_tsv(summary_path, required=["sample_id", "pseudogenized_count", "deleted_count", "uncertain_count"])
    if not rows:
        raise SchemaError(f"{summary_path}: no rows to plot")
    samples = [row["sample_id"] for row in rows]
    pseudogenized = [float(row["pseudogenized_count"]) for row in rows]
    deleted = [float(row["deleted_count"]) for row in rows]
    uncertain = [float(row["uncertain_count"]) for row in rows]
    x = list(range(len(samples)))
    figure, axis = plt.subplots(figsize=(max(8, len(samples) * 0.55), 6), constrained_layout=True)
    axis.bar(x, pseudogenized, label="Pseudogenized", color="#4d4d4d")
    axis.bar(x, deleted, bottom=pseudogenized, label="Deleted", color="#d9d9d9")
    axis.bar(x, uncertain, bottom=[a + b for a, b in zip(pseudogenized, deleted)], label="Uncertain", color="#e69f00")
    axis.set_ylabel("Reference genes")
    axis.set_xticks(x, samples, rotation=45, ha="right")
    if title:
        axis.set_title(title)
    axis.legend(frameon=False)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return output


def _numeric_chromosome(label: str) -> tuple[int, str]:
    digits = "".join(character for character in label if character.isdigit())
    return (int(digits), label) if digits else (10**9, label)


def plot_spatial_bubble(
    spatial_path: str | Path,
    output_path: str | Path,
    mode: str = "inter",
    rate_column: str = "loss_fragment_per_target_gene",
) -> Path:
    """Plot one dot per sample/chromosome (inter) or sample/bin/chromosome (intra)."""
    plt = _plt()
    required = ["sample_id", "target_chromosome", rate_column]
    if mode == "intra":
        required.append("bin")
    elif mode != "inter":
        raise SchemaError("mode must be inter or intra")
    rows = read_tsv(spatial_path, required=required)
    rows = [row for row in rows if row[rate_column] != ""]
    if not rows:
        raise SchemaError(f"{spatial_path}: no finite {rate_column} values to plot")
    samples = sorted({row["sample_id"] for row in rows})
    chromosomes = sorted({row["target_chromosome"] for row in rows}, key=_numeric_chromosome)
    sample_index = {sample: index for index, sample in enumerate(samples)}
    chromosome_index = {chromosome: index for index, chromosome in enumerate(chromosomes)}
    values = [float(row[rate_column]) for row in rows]
    if mode == "inter":
        x = [sample_index[row["sample_id"]] for row in rows]
        sizes = [80] * len(rows)
        x_label = "Sample"
    else:
        bins = sorted({int(row["bin"]) for row in rows})
        # A bin axis with small sample jitter prevents a categorical hue from
        # hiding same-bin points.
        offset = 0.7 / max(1, len(samples))
        x = [int(row["bin"]) + (sample_index[row["sample_id"]] - (len(samples) - 1) / 2) * offset for row in rows]
        maximum = max(values) or 1.0
        sizes = [30 + 470 * value / maximum for value in values]
        x_label = "Equal-width chromosome bin"
    figure, axis = plt.subplots(figsize=(max(9, len(samples) * 0.7), max(6, len(chromosomes) * 0.25)), constrained_layout=True)
    scatter = axis.scatter(
        x, [chromosome_index[row["target_chromosome"]] for row in rows], s=sizes, c=values,
        cmap="YlOrRd", edgecolor="black", linewidth=0.3, alpha=0.9,
    )
    colorbar = figure.colorbar(scatter, ax=axis, pad=0.02)
    colorbar.set_label("Pseudogene fragments per target gene")
    if mode == "inter":
        axis.set_xticks(range(len(samples)), samples, rotation=45, ha="right")
    else:
        axis.set_xticks(bins, [f"Bin {item}" for item in bins])
    axis.set_yticks(range(len(chromosomes)), chromosomes)
    axis.set_xlabel(x_label)
    axis.set_ylabel("Target chromosome")
    axis.grid(axis="y", linestyle="--", alpha=0.25)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
    return output
