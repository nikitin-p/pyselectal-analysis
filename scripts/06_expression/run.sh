#!/usr/bin/env bash
# Stage 11: expression dynamics of interesting TSS clusters across conditions.
#
# Prerequisites:
#   - Stage 8+9 must have run (results/tss/tss_matrix.tsv, tss_clustered.tsv)
#   - Stage 10 optional but recommended (results/cager/tagclusters.tsv.gz for IQW plot)
#
# Usage:
#   bash scripts/06_expression/run.sh [--force]

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

FIGURES_DIR="$(yaml_get results_figures "$PARAMS")"
FIGURES_DIR="${FIGURES_DIR:-results/figures}"
CAGER_DIR="$(yaml_get results_cager "$PARAMS")"
CAGER_DIR="${CAGER_DIR:-results/cager}"

log_info "Stage 11: expression dynamics"

$RSCRIPT "$SCRIPT_DIR/01_expression_dynamics.R" \
    --matrix      "$ROOT_DIR/results/tss/tss_matrix.tsv" \
    --clustered   "$ROOT_DIR/results/tss/tss_clustered.tsv" \
    --tagclusters "$ROOT_DIR/$CAGER_DIR/tagclusters.tsv.gz" \
    --samples     "$SAMPLES" \
    --outdir      "$ROOT_DIR/$FIGURES_DIR" \
    --params      "$PARAMS" \
    $FORCE_FLAG

log_info "Stage 11 complete."
log_info "  figures → $FIGURES_DIR/11_*.pdf"
