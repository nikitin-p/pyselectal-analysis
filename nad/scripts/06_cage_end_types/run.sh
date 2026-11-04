#!/usr/bin/env bash
# NAD branch step 6: CAGE end types at NAD captureSeq 5' ends.
#
# Runs on the Mac, where results/tss/typed_ctss.tsv.gz lives. Needs the two
# step-5 tables copied from vavilov into nad/results/05_deseq/:
#   scp -P 22 'pnikitin@85.89.112.24:~/pyselectal-analysis/nad/results/05_deseq/*_S_vs_N.tsv.gz' nad/results/05_deseq/
# Outputs go to nad/results/06_cage_end_types/ and
# nad/figures/06_end_types_nad_vs_control.pdf.
# Spec: docs/specs/2026-11-04-nad-step6-cage-end-types.md
#
# Usage:
#   bash nad/scripts/06_cage_end_types/run.sh [--params config/params.yaml] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
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
# package this step needs, so fall back to it rather than failing.
CONDA_ENV="$(yaml_get conda_env_r "$PARAMS")"
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    RSCRIPT=$(conda run -n "$CONDA_ENV" which Rscript)
    log_info "using Rscript from conda env: $CONDA_ENV"
else
    require_cmd Rscript
    RSCRIPT=$(command -v Rscript)
    log_warn "conda env '$CONDA_ENV' not found; using system Rscript"
fi

RESULTS_ROOT="$(yaml_get results_root "$PARAMS")"
FASTA="$(yaml_get genome_fasta "$PARAMS")"
ls "$NAD_DIR"/results/05_deseq/*_S_vs_N.tsv.gz > /dev/null 2>&1 \
    || die "no step-5 tables in nad/results/05_deseq/ (copy them from vavilov, see header)"
mkdir -p "$NAD_DIR/logs"
LOG="$NAD_DIR/logs/06_cage_end_types_$(date +%Y%m%d_%H%M%S).log"

log_info "NAD step 6: CAGE end types at NAD 5' ends"
cd "$ROOT_DIR"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/06_cage_end_types.R" \
    --ctss    "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
    --nad-dir "nad/results/05_deseq" \
    --fasta   "$FASTA" \
    --outdir  "nad/results/06_cage_end_types" \
    --figdir  "nad/figures" $FORCE_FLAG
