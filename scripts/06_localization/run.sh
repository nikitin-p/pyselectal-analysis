#!/usr/bin/env bash
# Stage 06: Genomic localization by 5' end type
#
# Compares where 1Sg, 2Sgg, 3Sggg, M reads localize across categories:
# enhancer / promoter / genic / intergenic.
#
# Usage:
#   bash scripts/06_localization/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2";  shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --force)   FORCE_FLAG="--force"; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"

# Genome and annotation come from the params file, so that switching datasets is
# a change of params file rather than a change of script. Relative paths are
# resolved against the project root; absolute paths are used as given.
resolve_path() {
    local v="$1"
    [[ -z "$v" ]] && { printf ''; return; }
    [[ "$v" != /* ]] && v="$ROOT_DIR/$v"
    printf '%s' "$v"
}
ANNOTATION_GTF="$(resolve_path "$(yaml_get annotation_gtf "$PARAMS")")"

CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)

log_info "Stage 06: Genomic localization by end type"

$RSCRIPT "$SCRIPT_DIR/01_localization.R" \
    --samples    "$SAMPLES" \
    --params     "$PARAMS" \
    --bam_sel    "$RESULTS_ROOT/bam_selected" \
    --bam_sub    "$RESULTS_ROOT/bam_subsampled" \
    --tss        "$RESULTS_ROOT/tss/tss_clustered.tsv" \
    --gtf        "$ANNOTATION_GTF" \
    --outdir     "$RESULTS_ROOT/localization" \
    --figdir     "$RESULTS_ROOT/figures" \
    $FORCE_FLAG

log_info "Stage 06 complete."
log_info "  tables  → results/localization/"
log_info "  figures → results/figures/10_*.pdf"
