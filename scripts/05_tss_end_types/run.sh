#!/usr/bin/env bash
# Stage 8: TSS × end-type matrix — per-cluster 1Sg/2Sgg/M/other fractions
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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2";  shift 2 ;;
        --samples) SAMPLES="$2"; shift 2 ;;
        --force)   FORCE_FLAG="--force"; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"

CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)

log_info "Stage 8: TSS end-type matrix"

log_info "  step 1: build typed CTSS"
$RSCRIPT "$SCRIPT_DIR/01_build_ctss.R" \
    --samples "$SAMPLES" \
    --outdir  "$RESULTS_ROOT/tss" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 2: cluster + annotate"
$RSCRIPT "$SCRIPT_DIR/02_cluster_annotate.R" \
    --ctss    "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
    --samples "$SAMPLES" \
    --outdir  "$RESULTS_ROOT/tss" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 3: plot"
$RSCRIPT "$SCRIPT_DIR/03_plot.R" \
    --matrix  "$RESULTS_ROOT/tss/tss_matrix.tsv" \
    --samples "$SAMPLES" \
    --outdir  "$RESULTS_ROOT/figures" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 4: cluster + annotate"
$RSCRIPT "$SCRIPT_DIR/04_cluster_annotate.R" \
    --matrix  "$RESULTS_ROOT/tss/tss_matrix.tsv" \
    --samples "$SAMPLES" \
    --outdir  "$RESULTS_ROOT/tss" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 5: plot clusters"
$RSCRIPT "$SCRIPT_DIR/05_plot_clusters.R" \
    --clustered  "$RESULTS_ROOT/tss/tss_clustered.tsv" \
    --hclust_rds "$RESULTS_ROOT/tss/tss_hclust.rds" \
    --outdir     "$RESULTS_ROOT/figures" \
    $FORCE_FLAG

# Relative paths in the params file are resolved against the project root.
ANNOTATION_BED="$(yaml_get annotation_bed "$PARAMS")"
[[ -n "$ANNOTATION_BED" && "$ANNOTATION_BED" != /* ]] && ANNOTATION_BED="$ROOT_DIR/$ANNOTATION_BED"

if [[ "$(yaml_get promoter_source "$PARAMS")" == "clusters" ]]; then
    log_info "  step 12: cluster-defined promoters"
    $RSCRIPT "$SCRIPT_DIR/12_promoter_clusters.R" \
        --ctss     "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
        --matrix   "$RESULTS_ROOT/tss/tss_matrix.tsv" \
        --genes    "$ANNOTATION_BED" \
        --outdir   "$RESULTS_ROOT/tss" \
        --trackdir "$RESULTS_ROOT/tracks" \
        --params   "$PARAMS" \
        $FORCE_FLAG

    log_info "  step 13: promoter window sweep and definition figure"
    $RSCRIPT "$SCRIPT_DIR/13_promoter_window_sweep.R" \
        --ctss     "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
        --matrix   "$RESULTS_ROOT/tss/tss_matrix.tsv" \
        --genes    "$ANNOTATION_BED" \
        --outdir   "$RESULTS_ROOT/tss" \
        --figdir   "$RESULTS_ROOT/figures" \
        --params   "$PARAMS" \
        $FORCE_FLAG
fi

log_info "Stage 8+9 complete."
log_info "  matrix    → results/tss/tss_matrix.tsv"
log_info "  clustered → results/tss/tss_clustered.tsv"
log_info "  figures   → results/figures/08_*.pdf  09_*.pdf"
