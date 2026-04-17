#!/usr/bin/env python3
"""Create a gene-activity matrix from a SnapATAC2 h5ad object."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Mapping, Optional

import scanpy as sc
import snapatac2 as snap

LOGGER = logging.getLogger(__name__)

GENOME_GFF_ENV: Mapping[str, str] = {
    "mm10": "MM10_GFF3",
    "hg38": "HG38_GFF3",
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create a SnapATAC2 gene-activity h5ad from a spatial ATAC h5ad."
    )
    parser.add_argument("-i", "--input_h5ad", required=True, help="Input ATAC h5ad file.")
    parser.add_argument("--species", default="mm10", choices=sorted(GENOME_GFF_ENV), help="Genome build.")
    parser.add_argument("--gff", default=None, help="Genome annotation GFF3 file. Overrides species env vars.")
    parser.add_argument("--min_cells", default=10, type=int, help="Minimum cells per gene.")
    parser.add_argument("--min_genes", default=200, type=int, help="Minimum genes per cell.")
    parser.add_argument("--output_h5ad", default=None, help="Output gene-activity h5ad file.")
    return parser.parse_args()


def resolve_gff(species: str, gff: Optional[str]) -> str:
    """Resolve the gene annotation file for a species.

    Args:
        species: Genome build, currently ``mm10`` or ``hg38``.
        gff: Optional explicit GFF3 path from the CLI.

    Returns:
        Path to a GFF3 annotation file.

    Raises:
        ValueError: If no annotation path is available for the species.
        FileNotFoundError: If the resolved annotation file is missing.
    """
    annotation = gff or os.environ.get(GENOME_GFF_ENV[species])
    if not annotation:
        raise ValueError(
            f"No GFF3 annotation configured for {species}. Set {GENOME_GFF_ENV[species]} "
            "or pass --gff."
        )

    annotation_path = Path(annotation)
    if not annotation_path.is_file():
        raise FileNotFoundError(f"GFF3 annotation file not found: {annotation_path}")
    return str(annotation_path)


def create_gene_activity(args: argparse.Namespace) -> Path:
    """Create and write a filtered gene-activity h5ad file.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Output h5ad path.
    """
    input_h5ad = Path(args.input_h5ad)
    if not input_h5ad.is_file():
        raise FileNotFoundError(f"Input h5ad file not found: {input_h5ad}")

    output_h5ad = Path(args.output_h5ad) if args.output_h5ad else input_h5ad.with_suffix(".gene_activity.h5ad")
    annotation = resolve_gff(args.species, args.gff)

    LOGGER.info("Reading ATAC h5ad: %s", input_h5ad)
    atac = snap.read(str(input_h5ad), backed=None)
    LOGGER.info("Creating gene-activity matrix with annotation: %s", annotation)
    gene_activity = snap.pp.make_gene_matrix(atac, annotation)

    sc.pp.filter_genes(gene_activity, min_cells=args.min_cells)
    sc.pp.filter_cells(gene_activity, min_genes=args.min_genes)

    LOGGER.info("Writing gene-activity h5ad: %s", output_h5ad)
    gene_activity.write(output_h5ad, compression="gzip")
    return output_h5ad


def main() -> None:
    """Run the gene-activity CLI."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    args = parse_args()
    create_gene_activity(args)


if __name__ == "__main__":
    main()
