#!/usr/bin/env bash
# Stage 17: TSS clusters dominated by a single 5' end type.
#
# Steps:
#   1. Re-classify reads by their literal 5' soft-clip sequence
#   2. Re-derive the Stage 8 clusters, assert identity, call dominance
#   3. Compare expression across dominance classes
#   4. Genomic context: templated artefact and G-masking controls
#
# Usage:
#   bash scripts/10_dominant_endtype/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
for arg in "$@"; do
    if [[ "$arg" == "--force" ]]; then
        FORCE_FLAG="--force"
    fi
done

# The r_cage environment may be absent; the system R library carries every
# package this stage needs, so fall back to it rather than failing.
CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)
    log_info "using Rscript from conda env: $CONDA_ENV"
else
    require_cmd Rscript
    RSCRIPT=$(command -v Rscript)
    log_warn "conda env '$CONDA_ENV' not found; using system Rscript"
fi

DOM_DIR="$ROOT_DIR/$(yaml_get results_dominance "$PARAMS")"

LOG="$ROOT_DIR/logs/17_dominance_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 17: single-type-dominated TSSs"

log_info "  step 1: fine-typed CTSS"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_build_fine_ctss.R" \
    --samples "$SAMPLES" \
    --outdir  "$DOM_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 2: clusters + dominance cascade"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_dominance.R" \
    --ctss    "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --matrix  "$ROOT_DIR/results/tss/tss_matrix.tsv" \
    --outdir  "$DOM_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 3: expression by class"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/03_expression.R" \
    --calls   "$DOM_DIR/dominance_calls.tsv" \
    --outdir  "$DOM_DIR" \
    --figdir  "$ROOT_DIR/results/figures" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 4: genomic context"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/04_context.R" \
    --calls   "$DOM_DIR/dominance_calls.tsv" \
    --matrix  "$ROOT_DIR/results/tss/tss_matrix.tsv" \
    --samples "$SAMPLES" \
    --outdir  "$DOM_DIR" \
    --figdir  "$ROOT_DIR/results/figures" $FORCE_FLAG

log_info "Stage 17 complete."
log_info "  calls    → results/dominance/dominance_calls.tsv"
log_info "  cascade  → results/dominance/filter_cascade.tsv"
log_info "  aggregate → results/dominance/aggregate_unusual.tsv"
log_info "  context  → results/dominance/templated_fraction.tsv"
log_info "  figures  → results/figures/18_dom_*.pdf"
log_info "  log      → $LOG"
