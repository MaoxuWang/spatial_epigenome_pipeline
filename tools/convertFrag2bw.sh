#!/usr/bin/env bash

# Convert a fragments.tsv.gz file to bedGraph and bigWig coverage tracks.
# The input is expected to contain at least chr, start, end in the first three columns.

set -euo pipefail

usage() {
    cat <<USAGE
Usage:
  bash tools/convertFrag2bw.sh <fragments.tsv.gz> <out_prefix> <chrom.sizes>

Environment:
  BEDTOOLS             Path to bedtools. Default: bedtools
  BEDGRAPH_TO_BIGWIG   Path to bedGraphToBigWig. Default: bedGraphToBigWig
USAGE
}

if [[ $# -ne 3 ]]; then
    usage >&2
    exit 1
fi

fragments_file="$1"
out_prefix="$2"
chrom_size="$3"
bedtools="${BEDTOOLS:-bedtools}"
bedgraph_to_bigwig="${BEDGRAPH_TO_BIGWIG:-bedGraphToBigWig}"

if [[ ! -s "$fragments_file" ]]; then
    echo "ERROR: fragments file does not exist or is empty: $fragments_file" >&2
    exit 1
fi
if [[ ! -s "$chrom_size" ]]; then
    echo "ERROR: chromosome sizes file does not exist or is empty: $chrom_size" >&2
    exit 1
fi

echo "[INFO] Writing bedGraph: ${out_prefix}.bedGraph" >&2
gzip -dc "$fragments_file" \
    | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2, $3}' \
    | "$bedtools" genomecov -bg -i stdin -g "$chrom_size" \
    > "${out_prefix}.bedGraph"

echo "[INFO] Sorting bedGraph: ${out_prefix}.sorted.bedGraph" >&2
sort -k1,1 -k2,2n "${out_prefix}.bedGraph" > "${out_prefix}.sorted.bedGraph"

echo "[INFO] Writing bigWig: ${out_prefix}.bw" >&2
"$bedgraph_to_bigwig" "${out_prefix}.sorted.bedGraph" "$chrom_size" "${out_prefix}.bw"

rm -f "${out_prefix}.bedGraph" "${out_prefix}.sorted.bedGraph"
echo "[INFO] Done: ${out_prefix}.bw" >&2
