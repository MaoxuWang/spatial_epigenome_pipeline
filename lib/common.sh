#!/usr/bin/env bash

# Shared shell helpers for the spatial epigenome pipelines.
# Pipeline entrypoints source this file before defining tool paths.

set -o pipefail

resolve_repo_root() {
    local source_file="$1"
    local source_dir
    source_dir="$(cd "$(dirname "$source_file")" && pwd)"
    cd "${source_dir}/.." && pwd
}

load_pipeline_config() {
    local repo_root="$1"
    local config_file="${SPATIAL_EPIGENOME_CONFIG:-${repo_root}/config/pipeline.env}"

    if [[ -f "$config_file" ]]; then
        # shellcheck source=/dev/null
        source "$config_file"
    fi
}

threads_or_default() {
    local default_threads="${1:-4}"
    local requested_threads="${SLURM_CPUS_PER_TASK:-$default_threads}"

    if [[ "$requested_threads" =~ ^[0-9]+$ ]] && [[ "$requested_threads" -gt 0 ]]; then
        echo "$requested_threads"
    else
        echo "$default_threads"
    fi
}

min_threads() {
    local requested="$1"
    local maximum="$2"

    if [[ "$requested" -gt "$maximum" ]]; then
        echo "$maximum"
    else
        echo "$requested"
    fi
}

require_file() {
    local path="$1"
    local label="$2"

    if [[ ! -s "$path" ]]; then
        echo "ERROR: ${label} does not exist or is empty: ${path}" >&2
        return 1
    fi
}

log_info() {
    echo "[INFO] $*" >&2
}

log_error() {
    echo "[ERROR] $*" >&2
}
