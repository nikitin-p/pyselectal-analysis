#!/usr/bin/env bash
# Stage 19: origins of the exotic 5' end types.
#
# "Exotic" means every fine 5'-end type other than 1Sg, so the unclipped class M
# is included. The stage asks where each exotic type's signal comes from, in two
# currencies at once (alignments and distinct CTSS positions), and how
# productive the positions it comes from are.
#
# Steps:
#   1. Re-run the Stage 17 dominance cascade with the CTSS position as the unit
#   2. Tally the exotics by origin and draw the composite figure
#
# This stage reads Stage 17's tables only. Run scripts/10_dominant_endtype/run.sh
# first; no BAM is opened here and nothing is re-typed.
#
# Usage:
#   bash scripts/12_exotic_origins/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
FORCE_FLAG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params) PARAMS="$2"; shift 2 ;;
        --force)  FORCE_FLAG="--force"; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
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
EXO_DIR="$ROOT_DIR/$(yaml_get results_exotic_origins "$PARAMS")"
FIG_DIR="$ROOT_DIR/$(yaml_get results_figures "$PARAMS")"

LOG="$ROOT_DIR/logs/19_exotic_origins_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 19: origins of the exotic 5' end types"

log_info "  step 1: position-level dominance cascade"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_position_dominance.R" \
    --ctss       "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --background "$DOM_DIR/background_shares.tsv"  \
    --outdir     "$EXO_DIR" \
    --params     "$PARAMS" $FORCE_FLAG

log_info "  step 2: tallies and composite figure"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_origins_figure.R" \
    --dominance      "$EXO_DIR/position_dominance.tsv.gz" \
    --ctss           "$DOM_DIR/fine_typed_ctss.tsv.gz"    \
    --cluster_counts "$DOM_DIR/cluster_type_counts.tsv"   \
    --outdir         "$EXO_DIR" \
    --figdir         "$FIG_DIR" \
    --params         "$PARAMS" $FORCE_FLAG

log_info "Stage 19 complete. Tables: $EXO_DIR  Figure: $FIG_DIR/20_exo_origins.pdf"
