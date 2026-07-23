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
[[ "${1:-}" == "--force" ]] && FORCE_FLAG="--force"

CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)

log_info "Stage 06: Genomic localization by end type"

$RSCRIPT "$SCRIPT_DIR/01_localization.R" \
    --samples    "$SAMPLES" \
    --params     "$PARAMS" \
    --bam_sel    "$ROOT_DIR/results/bam_selected" \
    --bam_sub    "$ROOT_DIR/results/bam_subsampled" \
    --tss        "$ROOT_DIR/results/tss/tss_clustered.tsv" \
    --gtf        "$ROOT_DIR/data/sacCer3.ensGene.gtf" \
    --outdir     "$ROOT_DIR/results/localization" \
    --figdir     "$ROOT_DIR/results/figures" \
    $FORCE_FLAG

log_info "Stage 06 complete."
log_info "  tables  → results/localization/"
log_info "  figures → results/figures/10_*.pdf"
