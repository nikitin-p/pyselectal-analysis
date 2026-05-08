#!/usr/bin/env bash
# Stage 10: CAGEr promoter clusters — normalize, distclu, IQW, annotation, plots.
#
# Steps:
#   1. Build CAGEexp from typed_ctss.tsv.gz → distclu + IQW → RDS + tagclusters TSV
#   2. Annotate consensus tag clusters against gene annotation BED
#   3. Diagnostic plots: correlation, reverse cumulatives, IQW, annotation enrichment
#
# Prerequisites:
#   - Stage 8 must have run (results/tss/typed_ctss.tsv.gz must exist)
#
# Usage:
#   bash scripts/06_cager/run.sh [--force]

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

CAGER_DIR="$(yaml_get results_cager "$PARAMS")"
CAGER_DIR="${CAGER_DIR:-results/cager}"
FIGURES_DIR="$(yaml_get results_figures "$PARAMS")"
FIGURES_DIR="${FIGURES_DIR:-results/figures}"

log_info "Stage 10: CAGEr promoter clusters"

log_info "  step 1: build CAGEexp + distclu + IQW"
$RSCRIPT "$SCRIPT_DIR/01_cageexp.R" \
    --ctss    "$ROOT_DIR/results/tss/typed_ctss.tsv.gz" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/$CAGER_DIR" \
    --params  "$PARAMS" \
    $FORCE_FLAG

log_info "  step 2: annotate tag clusters"
$RSCRIPT "$SCRIPT_DIR/02_annotate.R" \
    --tagclusters "$ROOT_DIR/$CAGER_DIR/tagclusters.tsv.gz" \
    --samples     "$SAMPLES" \
    --outdir      "$ROOT_DIR/$CAGER_DIR" \
    --params      "$PARAMS" \
    $FORCE_FLAG

log_info "  step 3: diagnostic plots"
$RSCRIPT "$SCRIPT_DIR/03_plot.R" \
    --rds          "$ROOT_DIR/$CAGER_DIR/cageexp.rds" \
    --annot        "$ROOT_DIR/$CAGER_DIR/annotation_enrichment.tsv" \
    --tagclusters  "$ROOT_DIR/$CAGER_DIR/tagclusters_annotated.tsv" \
    --outdir       "$ROOT_DIR/$FIGURES_DIR" \
    --params       "$PARAMS" \
    $FORCE_FLAG

log_info "Stage 10 complete."
log_info "  cageexp    → $CAGER_DIR/cageexp.rds"
log_info "  clusters   → $CAGER_DIR/tagclusters_annotated.tsv"
log_info "  figures    → $FIGURES_DIR/10_*.pdf"
