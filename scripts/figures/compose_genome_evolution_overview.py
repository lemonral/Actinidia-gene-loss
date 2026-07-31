#!/usr/bin/env python3
"""Combine four completed publication panels without rerunning analyses."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Sequence

from PIL import Image, ImageChops, ImageDraw, ImageFont


SCRIPT_VERSION = "1.0.0"
PANEL_LABELS = ("(a)", "(b)", "(c)", "(d)")


class CompositionError(RuntimeError):
    """Raised when source panels cannot be combined safely."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def trim_white_margin(image: Image.Image, *, padding: int = 16) -> Image.Image:
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "white")
    difference = ImageChops.difference(rgb, background).convert("L")
    difference = difference.point(lambda value: 255 if value > 3 else 0)
    bounds = difference.getbbox()
    if bounds is None:
        raise CompositionError("a source panel is blank")
    left, top, right, bottom = bounds
    left = max(0, left - padding)
    top = max(0, top - padding)
    right = min(rgb.width, right + padding)
    bottom = min(rgb.height, bottom + padding)
    return rgb.crop((left, top, right, bottom))


def publication_font(size: int, requested: Path | None) -> ImageFont.FreeTypeFont:
    candidates = [
        requested,
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        Path("/System/Library/Fonts/Supplemental/Arial.ttf"),
        Path("/Library/Fonts/Arial.ttf"),
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return ImageFont.truetype(str(candidate), size=size)
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size=size)
    except OSError as error:
        raise CompositionError(
            "no publication font was found; pass --font with a TrueType font"
        ) from error


def fit_panel(image: Image.Image, width: int, height: int) -> Image.Image:
    scale = min(width / image.width, height / image.height)
    if scale <= 0:
        raise CompositionError("invalid panel dimensions")
    target = (
        max(1, round(image.width * scale)),
        max(1, round(image.height * scale)),
    )
    return image.resize(target, Image.Resampling.LANCZOS)


def compose(
    *,
    sources: Sequence[Path],
    output_dir: Path,
    basename: str,
    canvas_size: int,
    dpi: int,
    font_path: Path | None,
) -> dict[str, object]:
    if len(sources) != 4:
        raise CompositionError("exactly four source panels are required")
    for source in sources:
        if not source.is_file():
            raise CompositionError(f"missing source panel: {source}")
    if canvas_size < 2400:
        raise CompositionError("--canvas-size must be at least 2400 pixels")
    if dpi < 150:
        raise CompositionError("--dpi must be at least 150")

    output_dir.mkdir(parents=True, exist_ok=True)
    outer = round(canvas_size * 0.018)
    gap = round(canvas_size * 0.014)
    label_band = round(canvas_size * 0.022)
    cell_width = (canvas_size - 2 * outer - gap) // 2
    cell_height = (canvas_size - 2 * outer - gap) // 2
    image_width = cell_width
    image_height = cell_height - label_band
    label_font = publication_font(round(canvas_size * 0.016), font_path)

    canvas = Image.new("RGB", (canvas_size, canvas_size), "white")
    draw = ImageDraw.Draw(canvas)
    source_rows: list[dict[str, object]] = []

    for index, source in enumerate(sources):
        row, column = divmod(index, 2)
        x0 = outer + column * (cell_width + gap)
        y0 = outer + row * (cell_height + gap)
        panel = trim_white_margin(Image.open(source))
        panel = fit_panel(panel, image_width, image_height)
        paste_x = x0 + (image_width - panel.width) // 2
        paste_y = y0 + label_band + (image_height - panel.height) // 2
        canvas.paste(panel, (paste_x, paste_y))
        draw.text((x0, y0), PANEL_LABELS[index], fill="black", font=label_font)
        source_rows.append(
            {
                "panel": PANEL_LABELS[index],
                "basename": source.name,
                "bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )

    png_path = output_dir / f"{basename}.png"
    pdf_path = output_dir / f"{basename}.pdf"
    caption_path = output_dir / f"{basename}.caption.txt"
    manifest_path = output_dir / f"{basename}.manifest.json"

    canvas.save(png_path, dpi=(dpi, dpi), optimize=True)
    canvas.save(pdf_path, "PDF", resolution=dpi, quality=95)
    caption = (
        "(a) Chromosome-scale features and intragenomic syntenic connections "
        "in Clematoclethra scandens. (b) OrthoFinder gene-category composition "
        "for the frozen 17-taxon dataset. (c) Kernel-density distributions of "
        "synonymous substitutions per synonymous site (Ks). (d) The "
        "TimeTree-secondary-calibrated species tree with CAFE5 Base-model "
        "branch expansions and contractions.\n"
    )
    caption_path.write_text(caption, encoding="utf-8")

    outputs = []
    for path, role in (
        (png_path, "figure_png"),
        (pdf_path, "figure_pdf"),
        (caption_path, "caption_text"),
    ):
        outputs.append(
            {
                "basename": path.name,
                "role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "script_version": SCRIPT_VERSION,
        "layout": "two_by_two",
        "canvas_pixels": [canvas_size, canvas_size],
        "dpi": dpi,
        "inputs": source_rows,
        "outputs": outputs,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--circos", required=True, type=Path)
    parser.add_argument("--orthofinder", required=True, type=Path)
    parser.add_argument("--ks", required=True, type=Path)
    parser.add_argument("--phylogeny", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="genome_evolution_overview")
    parser.add_argument("--canvas-size", type=int, default=4800)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font", type=Path)
    args = parser.parse_args(argv)
    try:
        manifest = compose(
            sources=(args.circos, args.orthofinder, args.ks, args.phylogeny),
            output_dir=args.output_dir,
            basename=args.basename,
            canvas_size=args.canvas_size,
            dpi=args.dpi,
            font_path=args.font,
        )
    except (OSError, CompositionError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
