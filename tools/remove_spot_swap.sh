#!/usr/bin/env bash

# Remove likely spot-swap / bleed-over fragments with a two-pass hybrid filter.
# Expected input format: chr, start, end, barcode, count, strand.
# Barcode names are expected to contain X/Y bins as fields 3 and 4 after splitting
# by underscores, for example: s_016um_00365_00044.

set -euo pipefail

usage() {
    cat <<USAGE
Usage:
  bash tools/remove_spot_swap.sh <input.fragments.tsv.gz> <output.cleaned.fragments.tsv.gz> [resolution_um] [local_radius_um] [local_threshold]

Defaults:
  resolution_um      16
  local_radius_um    64
  local_threshold    0.1

Notes:
  - The script requires awk with arrays-of-arrays support, bgzip, tabix, and sort.
  - Fragments within local_radius_um are kept if count >= dominant_count * local_threshold.
  - Fragments outside local_radius_um are kept only at the dominant barcode coordinate.
USAGE
}

if [[ $# -lt 2 || $# -gt 5 ]]; then
    usage >&2
    exit 1
fi

input_file="$1"
output_file="$2"
resolution_um="${3:-16}"
local_radius_um="${4:-64}"
local_threshold="${5:-0.1}"

if [[ ! -s "$input_file" ]]; then
    echo "ERROR: input fragments file does not exist or is empty: $input_file" >&2
    exit 1
fi

if [[ "$input_file" == *.gz ]]; then
    reader=(gzip -dc)
else
    reader=(cat)
fi

local_radius_bins_sq=$(awk -v radius="$local_radius_um" -v resolution="$resolution_um" 'BEGIN{r = radius / resolution; print r * r}')
tmp_map="${output_file}.dominant_map.$$.$RANDOM.tmp"
trap 'rm -f "$tmp_map"' EXIT

echo "[INFO] Pass 1/2: building dominant-barcode map" >&2
"${reader[@]}" "$input_file" | awk '
BEGIN {
    FS = "\t";
    OFS = "\t";
    format_checked = 0;
}
!/^#/ {
    if (format_checked == 0) {
        if (NF != 6) {
            printf("ERROR: expected 6 columns, found %d at line %d\n", NF, NR) > "/dev/stderr";
            exit 1;
        }
        format_checked = 1;
    }

    count = $5;
    barcode = $4;
    fragment_key = $1":"$2":"$3":"$6;

    if (!(barcode in barcode_to_xy)) {
        split(barcode, parts, "_");
        if (length(parts[3]) == 0 || length(parts[4]) == 0) {
            printf("ERROR: cannot parse X/Y coordinates from barcode %s\n", barcode) > "/dev/stderr";
            exit 1;
        }
        barcode_to_xy[barcode] = parts[3]"\t"parts[4];
    }

    fragment_counts_per_barcode[fragment_key][barcode] += count;
    current_count = fragment_counts_per_barcode[fragment_key][barcode];
    if (current_count > fragment_max_count[fragment_key]) {
        fragment_max_count[fragment_key] = current_count;
        dominant_barcode_for_fragment[fragment_key] = barcode;
    }
}
END {
    for (fragment_key in dominant_barcode_for_fragment) {
        dominant_barcode = dominant_barcode_for_fragment[fragment_key];
        print fragment_key, barcode_to_xy[dominant_barcode], fragment_max_count[fragment_key];
    }
}
' > "$tmp_map"

if [[ ! -s "$tmp_map" ]]; then
    echo "ERROR: dominant-barcode map is empty: $tmp_map" >&2
    exit 1
fi

echo "[INFO] Pass 2/2: applying hybrid local/remote filter" >&2
"${reader[@]}" "$input_file" | awk \
    -v radius_sq="$local_radius_bins_sq" \
    -v threshold="$local_threshold" \
    'BEGIN { FS = "\t"; OFS = "\t"; }
NR == FNR {
    dominant_x[$1] = $2;
    dominant_y[$1] = $3;
    dominant_count[$1] = $4;
    next;
}
/^#/ { next; }
{
    count = $5;
    barcode = $4;
    fragment_key = $1":"$2":"$3":"$6;

    split(barcode, parts, "_");
    current_x = parts[3];
    current_y = parts[4];
    dom_x = dominant_x[fragment_key];
    dom_y = dominant_y[fragment_key];
    dom_count = dominant_count[fragment_key];

    if (dom_x == "" || dom_y == "") {
        print $0;
        next;
    }

    dist_sq = (current_x - dom_x)^2 + (current_y - dom_y)^2;
    cutoff_count = dom_count * threshold;

    if (dist_sq <= radius_sq) {
        if (count >= cutoff_count) {
            print $0;
        }
    } else if (current_x == dom_x && current_y == dom_y) {
        print $0;
    }
}
' "$tmp_map" - \
    | sort -k1,1V -k2,2n -T . \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' \
    | bgzip > "$output_file"

tabix -p bed "$output_file"
echo "[INFO] Cleaned fragments written to $output_file" >&2
