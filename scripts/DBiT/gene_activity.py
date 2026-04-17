#!/usr/bin/env python3
"""Create DBiT-compatible gene activity outputs from h5ad or BAM input."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Mapping, Optional

import pandas as pd
import snapatac2 as snap

LOGGER = logging.getLogger(__name__)

GENOME_GFF_ENV: Mapping[str, str] = {
    "mm10": "MM10_GFF3",
    "hg38": "HG38_GFF3",
}
GENOME_SIZES = {
    "mm10": snap.genome.mm10,
    "hg38": snap.genome.hg38,
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create gene-activity h5ad/TSV files and spatial barcode mapping for DBiT data."
    )
    parser.add_argument("-i", "--input_h5ad", required=True, help="Input h5ad or BAM file.")
    parser.add_argument("--species", default="mm10", choices=sorted(GENOME_GFF_ENV), help="Genome build.")
    parser.add_argument("--gff", default=None, help="Genome annotation GFF3 file. Overrides species env vars.")
    parser.add_argument("--barcode_file", required=True, help="TSV with barcode, x, y columns.")
    parser.add_argument("--outdir", required=True, help="Output directory used by the calling pipeline.")
    parser.add_argument("--is_paired", action="store_true", help="Treat BAM as paired-end when making fragments.")
    parser.add_argument("--barcode_tag", default="CB", help="BAM tag containing the spatial barcode.")
    parser.add_argument("--output_prefix", default=None, help="Output prefix without suffix.")
    return parser.parse_args()


def resolve_gff(species: str, gff: Optional[str]) -> str:
    """Resolve and validate the GFF3 annotation file.

    Args:
        species: Genome build.
        gff: Optional explicit GFF3 path.

    Returns:
        Existing GFF3 path.

    Raises:
        ValueError: If the annotation is not configured.
        FileNotFoundError: If the annotation does not exist.
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


def read_spatial_positions(barcode_file: str) -> pd.DataFrame:
    """Read DBiT spatial barcode positions.

    Args:
        barcode_file: TSV with columns barcode, x, y and no header.

    Returns:
        Data frame indexed as ``<x>x<y>`` with a ``barcode`` column.

    Raises:
        ValueError: If the barcode file does not contain at least three columns.
    """
    barcodes = pd.read_csv(barcode_file, sep="\t", lineterminator="\n", header=None)
    if barcodes.shape[1] < 3:
        raise ValueError(f"Barcode file must have at least three columns: {barcode_file}")
    barcodes = barcodes.iloc[:, :3]
    barcodes.columns = ["barcode", "x", "y"]
    barcodes.index = barcodes["x"].astype(str) + "x" + barcodes["y"].astype(str)
    return barcodes.drop(columns=["x", "y"])


def read_or_make_atac_data(input_path: str, species: str, is_paired: bool, barcode_tag: str):
    """Read h5ad input or create a SnapATAC2 object from a BAM file.

    Args:
        input_path: h5ad or BAM path.
        species: Genome build.
        is_paired: Whether BAM input is paired-end.
        barcode_tag: BAM tag with spatial barcode.

    Returns:
        SnapATAC2/AnnData object for gene matrix creation.
    """
    if input_path.endswith(".h5ad"):
        LOGGER.info("Reading h5ad input: %s", input_path)
        return snap.read(input_path, backed=None)

    input_bam = Path(input_path)
    if not input_bam.is_file():
        raise FileNotFoundError(f"Input BAM/h5ad not found: {input_bam}")

    fragments = str(input_bam).replace(".bam", ".fragments.tsv.gz")
    LOGGER.info("Creating fragments from BAM: %s", input_bam)
    snap.pp.make_fragment_file(
        str(input_bam),
        fragments,
        is_paired=is_paired,
        barcode_tag=barcode_tag,
    )
    return snap.pp.import_data(
        fragments,
        chrom_sizes=GENOME_SIZES[species],
        sorted_by_barcode=False,
        min_num_fragments=0,
        n_jobs=-1,
    )


def write_outputs(args: argparse.Namespace) -> tuple[Path, Path]:
    """Create gene activity and write h5ad, expression TSV, and spatial barcodes.

    Args:
        args: Parsed command-line arguments.

    Returns:
        Tuple of output h5ad and expression TSV paths.
    """
    annotation = resolve_gff(args.species, args.gff)
    atac = read_or_make_atac_data(args.input_h5ad, args.species, args.is_paired, args.barcode_tag)
    LOGGER.info("Creating gene-activity matrix with annotation: %s", annotation)
    gene_activity = snap.pp.make_gene_matrix(atac, annotation)

    if args.output_prefix:
        output_prefix = Path(args.output_prefix)
    else:
        output_prefix = Path(args.input_h5ad.replace(".h5ad", ""))
    output_h5ad = output_prefix.with_suffix(".gene_activity.h5ad")
    output_tsv = output_prefix.with_suffix(".exp.tsv")
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    barcodes = read_spatial_positions(args.barcode_file)
    barcode_to_position = {row["barcode"]: position for position, row in barcodes.iterrows()}
    gene_activity.obs_names = [barcode_to_position[barcode.replace("-1", "")] for barcode in gene_activity.obs_names]

    spatial_dir = Path(args.outdir) / "spatial"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    barcodes.assign(index=barcodes.index).to_csv(spatial_dir / "spatial_barcodes.txt", index=False, sep="\t")

    LOGGER.info("Writing h5ad: %s", output_h5ad)
    gene_activity.write(output_h5ad, compression="gzip")

    matrix = gene_activity.X.toarray() if hasattr(gene_activity.X, "toarray") else gene_activity.X
    LOGGER.info("Writing expression TSV: %s", output_tsv)
    pd.DataFrame(matrix, index=gene_activity.obs_names, columns=gene_activity.var_names).to_csv(output_tsv, sep="\t")
    return output_h5ad, output_tsv


def main() -> None:
    """Run the DBiT gene-activity CLI."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    write_outputs(parse_args())


if __name__ == "__main__":
    main()
