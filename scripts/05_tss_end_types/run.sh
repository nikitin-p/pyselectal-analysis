#!/usr/bin/env bash
# Stage 8: TSS × end-type matrix — per-cluster 1Sg/2Sg/M/other fractions
# annotated with initiator dinucleotide.
#
# Steps:
#   1. Classify reads by 5'-end type and build typed CTSS table
#   2. distclu clustering (maxDist=20, CAGEr default) + initiator annotation
#   3. Figures: heatmap, initiator boxplot, stacked bar
#
# Usage:
#   bash scripts/05_tss_end_types/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
[[ "${1:-}" == "--force" ]] && FORCE_FLAG="--force"

CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT="conda run -n $CONDA_ENV Rscript"

log_info "Stage 8: TSS end-type matrix"

log_info "  step 1: build typed CTSS"
$RSCRIPT "$SCRIPT_DIR/01_build_ctss.R" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/results/tss" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 2: cluster + annotate"
$RSCRIPT "$SCRIPT_DIR/02_cluster_annotate.R" \
    --ctss    "$ROOT_DIR/results/tss/typed_ctss.tsv.gz" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/results/tss" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 3: plot"
$RSCRIPT "$SCRIPT_DIR/03_plot.R" \
    --matrix  "$ROOT_DIR/results/tss/tss_matrix.tsv" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/results/figures" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "Stage 8 complete."
log_info "  matrix  → results/tss/tss_matrix.tsv"
log_info "  figures → results/figures/08_*.pdf"
