#!/usr/bin/env bash

# Estimate spot-swap / bleed-over by counting fragments observed in more than one barcode.
# Expected fragments format: chr, start, end, barcode, count, strand.
# The implementation uses gawk arrays-of-arrays and can be memory intensive for very
# large fragment files because it indexes fragment -> barcode membership.

set -euo pipefail

usage() {
    cat <<USAGE
Usage:
  bash tools/cal_swap_score.sh <fragments.tsv[.gz]>
USAGE
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

fragments_file="$1"
if [[ ! -s "$fragments_file" ]]; then
    echo "ERROR: fragments file does not exist or is empty: $fragments_file" >&2
    exit 1
fi

if [[ "$fragments_file" == *.gz ]]; then
    reader=(gzip -dc)
else
    reader=(cat)
fi

echo "[INFO] Estimating spot-swap rate for $fragments_file" >&2

"${reader[@]}" "$fragments_file" | awk '
BEGIN {
    FS = "\t";
    OFS = "\t";
    total_reads = 0;
    total_swapped_reads = 0;
}
!/^#/ {
    if (NF < 6) {
        printf("ERROR: expected at least 6 columns, found %d at line %d\n", NF, NR) > "/dev/stderr";
        exit 1;
    }

    count = $5;
    fragment_key = $1":"$2":"$3":"$6;
    barcode = $4;

    total_reads += count;
    fragment_total_counts[fragment_key] += count;
    barcode_seen_for_fragment[fragment_key][barcode] = 1;
}
END {
    for (fragment_key in fragment_total_counts) {
        n_barcodes = 0;
        for (barcode in barcode_seen_for_fragment[fragment_key]) {
            n_barcodes++;
        }
        if (n_barcodes > 1) {
            total_swapped_reads += fragment_total_counts[fragment_key];
        }
    }

    if (total_reads > 0) {
        percentage = (total_swapped_reads / total_reads) * 100;
    } else {
        percentage = 0;
    }

    print "metric\tvalue";
    printf "total_reads\t%d\n", total_reads;
    printf "swapped_reads\t%d\n", total_swapped_reads;
    printf "spot_swap_percent\t%.4f\n", percentage;
}
'
