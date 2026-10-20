#!/usr/bin/env bash
# Stage 18: A-8 saturation of <5' end type, dinucleotide> pairs.
#
# Steps:
#   1. Build one compact window row per CTSS with the dinucleotide, the A-8
#      indicator, the two background A tallies and a distclu cluster label
#   2. Form the 88 testable pairs, test A-8 against each pair's own background,
#      and bootstrap the effect size over clusters
#   3. Render the bubble plots, the profile and the sensitivity comparison
#   4. Repeat the test with the dinucleotides pooled into YR / RY / RR / YY and
#      render the class-level bubble plots
#
# Usage:
#   bash scripts/11_a8_saturation/run.sh [--force]

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

A8_DIR="$ROOT_DIR/$(yaml_get results_a8 "$PARAMS")"
FIG_DIR="$ROOT_DIR/$(yaml_get results_figures "$PARAMS")"

LOG="$ROOT_DIR/logs/18_a8_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 18: A-8 saturation of <5' end type, dinucleotide> pairs"

log_info "  step 1: windows"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_build_windows.R" \
    --ctss    "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
    --samples "$SAMPLES" \
    --matrix  "$RESULTS_ROOT/tss/tss_matrix.tsv" \
    --outdir  "$A8_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 2: pair tests and bootstrap"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_a8_test.R" \
    --windows "$A8_DIR/pair_windows.tsv.gz" \
    --outdir  "$A8_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 3: figures"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/03_plot.R" \
    --tests   "$A8_DIR/a8_pair_tests.tsv" \
    --profile "$A8_DIR/a8_position_profile.tsv" \
    --figdir  "$FIG_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "  step 4: dinucleotide-class tests and figures"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/04_class_test.R" \
    --windows "$A8_DIR/pair_windows.tsv.gz" \
    --outdir  "$A8_DIR" \
    --figdir  "$FIG_DIR" \
    --params  "$PARAMS" $FORCE_FLAG

log_info "Stage 18 complete."
log_info "  tests   → results/a8/a8_pair_tests.tsv"
log_info "  classes → results/a8/a8_class_tests.tsv"
log_info "  sizes   → results/a8/a8_pair_sizes.tsv"
log_info "  profile → results/a8/a8_position_profile.tsv"
log_info "  figures → results/figures/19_a8_*.pdf"
log_info "  log     → $LOG"
