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

BAM_DIR="$ROOT_DIR/results/bam_subsampled"
SAMPLES="$ROOT_DIR/config/samples.tsv"
OUTDIR="$ROOT_DIR/results/dinuc"
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

log "Stage 6: dinucleotide proportions"
log "  bam_dir : $BAM_DIR"
log "  samples : $SAMPLES"
log "  outdir  : $OUTDIR"

mkdir -p "$OUTDIR"

Rscript "$SCRIPT_DIR/01_dinuc_proportion.R" \
    --bam_dir  "$BAM_DIR"   \
    --samples  "$SAMPLES"   \
    --outdir   "$OUTDIR"    \
    --params   "$PARAMS"    \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"

log "Stage 6 complete. Outputs in $OUTDIR"
