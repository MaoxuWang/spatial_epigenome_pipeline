#!/usr/bin/env python3
"""Split a RongFan spatial ATAC R1 read into sequence and barcode FASTQs."""

from __future__ import annotations

import argparse
import gzip
import logging
from pathlib import Path
from typing import TextIO

from Bio.SeqIO.QualityIO import FastqGeneralIterator

LOGGER = logging.getLogger(__name__)

SEQ_START = 99
BC2_SLICE = slice(0, 8)
BC1_SLICE = slice(38, 46)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create Cell Ranger ATAC compatible R1/R2 FASTQs from RongFan spatial ATAC R1."
    )
    parser.add_argument("-i", "--input", required=True, help="Input R1 FASTQ.GZ file.")
    parser.add_argument("-o1", "--output_R1", required=True, help="Output sequence FASTQ file.")
    parser.add_argument("-o2", "--output_R2", required=True, help="Output barcode FASTQ file.")
    return parser.parse_args()


def write_fastq_record(handle: TextIO, title: str, sequence: str, quality: str) -> None:
    """Write one FASTQ record.

    Args:
        handle: Writable text handle.
        title: FASTQ title without leading ``@``.
        sequence: Read sequence.
        quality: ASCII quality string.
    """
    handle.write(f"@{title}\n{sequence}\n+\n{quality}\n")


def process_fastq(input_fastq: Path, output_r1: Path, output_r2: Path) -> int:
    """Split sequence and barcode records from the input FASTQ.

    Args:
        input_fastq: Input gzipped FASTQ path.
        output_r1: Output sequence FASTQ path.
        output_r2: Output barcode FASTQ path.

    Returns:
        Number of processed reads.

    Raises:
        FileNotFoundError: If the input FASTQ is missing.
    """
    if not input_fastq.is_file():
        raise FileNotFoundError(f"Input FASTQ not found: {input_fastq}")

    output_r1.parent.mkdir(parents=True, exist_ok=True)
    output_r2.parent.mkdir(parents=True, exist_ok=True)

    n_reads = 0
    with gzip.open(input_fastq, "rt") as in_handle, output_r1.open("w") as out_r1, output_r2.open("w") as out_r2:
        for title, sequence, quality in FastqGeneralIterator(in_handle):
            new_sequence = sequence[SEQ_START:]
            new_quality = quality[SEQ_START:]
            barcode = sequence[BC2_SLICE] + sequence[BC1_SLICE]
            barcode_quality = quality[BC2_SLICE] + quality[BC1_SLICE]
            write_fastq_record(out_r1, title, new_sequence, new_quality)
            write_fastq_record(out_r2, title, barcode, barcode_quality)
            n_reads += 1

    return n_reads


def main() -> None:
    """Run the FASTQ splitting CLI."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    args = parse_args()
    n_reads = process_fastq(Path(args.input), Path(args.output_R1), Path(args.output_R2))
    LOGGER.info("Processed %d reads", n_reads)


if __name__ == "__main__":
    main()
