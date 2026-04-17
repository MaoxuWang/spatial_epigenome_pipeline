#!/usr/bin/env bash

# Split LIANTI/ME-like reads by detecting the Tn5/ME adaptor at the start of R1.
# This helper is intentionally small because the main pipelines call cutadapt directly
# for most production runs.

set -euo pipefail

usage() {
    cat <<USAGE
Usage:
  bash tools/cal_ME.sh <read1.fastq.gz> <read2.fastq.gz> <sample_id> <outdir>

Environment:
  CUTADAPT              Path to cutadapt executable. Default: cutadapt
  SLURM_CPUS_PER_TASK   Number of cutadapt threads. Default: 4
USAGE
}

if [[ $# -ne 4 ]]; then
    usage >&2
    exit 1
fi

read1="$1"
read2="$2"
sample_id="$3"
outdir="$4"
threads="${SLURM_CPUS_PER_TASK:-4}"
cutadapt="${CUTADAPT:-cutadapt}"

if [[ ! -s "$read1" ]]; then
    echo "ERROR: read1 does not exist or is empty: $read1" >&2
    exit 1
fi
if [[ ! -s "$read2" ]]; then
    echo "ERROR: read2 does not exist or is empty: $read2" >&2
    exit 1
fi

mkdir -p "${outdir}/trim"

"$cutadapt" \
    -j "$threads" \
    -g '^NAGATGTGTATAAGAGACAG' \
    -g '^NNAGATGTGTATAAGAGACAG' \
    -g '^NNNAGATGTGTATAAGAGACAG' \
    -g '^NNNNAGATGTGTATAAGAGACAG' \
    --no-indels \
    -e 2 \
    -o "${outdir}/trim/${sample_id}_ME_R1.fastq.gz" \
    -p "${outdir}/trim/${sample_id}_ME_R2.fastq.gz" \
    --untrimmed-output "${outdir}/trim/${sample_id}_woME_R1.fastq.gz" \
    --untrimmed-paired-output "${outdir}/trim/${sample_id}_woME_R2.fastq.gz" \
    "$read1" "$read2"
