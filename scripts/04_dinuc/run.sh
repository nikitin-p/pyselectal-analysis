#!/usr/bin/env bash
# Stage 6: initiator dinucleotide proportion per sample.
#
# Usage:
#   bash scripts/04_dinuc/run.sh [--bam_dir DIR] [--samples FILE] \
#                                [--outdir DIR] [--params FILE] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

BAM_DIR=""
SAMPLES="$ROOT_DIR/config/samples.tsv"
OUTDIR=""
PARAMS="$ROOT_DIR/config/params.yaml"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam_dir) BAM_DIR="$2";  shift 2 ;;
        --samples) SAMPLES="$2";  shift 2 ;;
        --outdir)  OUTDIR="$2";   shift 2 ;;
        --params)  PARAMS="$2";   shift 2 ;;
        --force)   EXTRA_ARGS+=(--force); shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# Defaults are resolved after argument parsing, so that they follow whichever
# params file was selected rather than the yeast tree.
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"
[[ -z "$BAM_DIR" ]] && BAM_DIR="$RESULTS_ROOT/bam_subsampled"
[[ -z "$OUTDIR"  ]] && OUTDIR="$RESULTS_ROOT/dinuc"

CONDA_R="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_R" which Rscript)

log_info "Stage 6: dinucleotide proportions"
log_info "  bam_dir : $BAM_DIR"
log_info "  samples : $SAMPLES"
log_info "  outdir  : $OUTDIR"

mkdir -p "$OUTDIR"

$RSCRIPT "$SCRIPT_DIR/01_dinuc_proportion.R" \
    --bam_dir  "$BAM_DIR"   \
    --samples  "$SAMPLES"   \
    --outdir   "$OUTDIR"    \
    --params   "$PARAMS"    \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"

log_info "Stage 6 complete. Outputs in $OUTDIR"
