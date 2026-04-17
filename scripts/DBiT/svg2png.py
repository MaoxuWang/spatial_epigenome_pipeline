#!/usr/bin/env python3
"""Convert an SVG file to PNG."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import cairosvg

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Convert one SVG file to PNG.")
    parser.add_argument("-i", "--input_svg", required=True, help="Input SVG file.")
    parser.add_argument("-o", "--output_png", required=True, help="Output PNG file.")
    return parser.parse_args()


def convert_svg_to_png(svg_file: Path, png_file: Path) -> None:
    """Convert SVG to PNG.

    Args:
        svg_file: Existing input SVG path.
        png_file: Output PNG path.

    Raises:
        FileNotFoundError: If the SVG input is missing.
    """
    if not svg_file.is_file():
        raise FileNotFoundError(f"Input SVG not found: {svg_file}")
    png_file.parent.mkdir(parents=True, exist_ok=True)
    cairosvg.svg2png(url=str(svg_file), write_to=str(png_file))
    LOGGER.info("Wrote PNG: %s", png_file)


def main() -> None:
    """Run the converter CLI."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    args = parse_args()
    convert_svg_to_png(Path(args.input_svg), Path(args.output_png))


if __name__ == "__main__":
    main()
