#!/usr/bin/env bash
# NAD branch step 8: NAD captureSeq signal and YAAG at CAGE single-A-clip positions.
#
# Runs on the Mac, where results/tss/typed_ctss.tsv.gz lives. Needs from vavilov:
#   scp -P 22 'pnikitin@85.89.112.24:~/pyselectal-analysis/nad/results/04_ctss/matrix_reads.tsv.gz' nad/results/04_ctss/
#   scp -P 22 'pnikitin@85.89.112.24:~/pyselectal-analysis/nad/results/05_deseq/*_S_vs_N.tsv.gz' nad/results/05_deseq/
# and, for the size-factor check, nad/results/07_chrM_sites/sanity.tsv (step 7).
# Outputs: nad/results/08_1Sa_positions/*.tsv, nad/figures/08_1Sa_positions.pdf.
# Spec: docs/specs/2026-11-04-nad-step8-1Sa-positions.md
#
# --self-test runs the two synthetic checks instead of the real analysis, each
# into its own folder so that the real tables are never overwritten:
#   selftest_null/   NAD values shuffled within strata (expect no excess)
#   selftest_plant/  30 % of clip-rich positions planted as NAD-detected (expect p ~ 0.001)
#
# Usage:
#   bash nad/scripts/08_1Sa_positions/run.sh [--params config/params.yaml] [--force] [--self-test]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
FORCE_FLAG=""; SELF=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)    PARAMS="$2"; shift 2 ;;
        --force)     FORCE_FLAG="--force"; shift ;;
        --self-test) SELF=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# Step 8 loads Biostrings and patchwork, which the system R library carries.
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
CTSS="$RESULTS_ROOT/tss/typed_ctss.tsv.gz"
MATRIX="nad/results/04_ctss/matrix_reads.tsv.gz"
cd "$ROOT_DIR"
[[ -s "$CTSS" ]]   || die "missing $CTSS (run Stage 8 of the main pipeline)"
[[ -s "$MATRIX" ]] || die "missing $MATRIX (copy it from vavilov, see header)"
ls nad/results/05_deseq/*_S_vs_N.tsv.gz > /dev/null 2>&1 \
    || die "no step-5 tables in nad/results/05_deseq/ (copy them from vavilov, see header)"
[[ -s nad/results/07_chrM_sites/sanity.tsv ]] \
    || log_warn "nad/results/07_chrM_sites/sanity.tsv missing: the size-factor check against step 7 is skipped"
mkdir -p "$NAD_DIR/logs"

run_step8() {   # run_step8 OUTDIR FIGDIR [extra args]
    local outdir="$1" figdir="$2"; shift 2
    log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/08_1Sa_positions.R" \
        --ctss "$CTSS" --matrix "$MATRIX" \
        --manifest "nad/scripts/00_download/priority1_runs.tsv" \
        --step5-dir "nad/results/05_deseq" --sanity "nad/results/07_chrM_sites/sanity.tsv" \
        --fasta "$FASTA" --outdir "$outdir" --figdir "$figdir" $FORCE_FLAG "$@"
}

if [[ $SELF -eq 1 ]]; then
    LOG="$NAD_DIR/logs/08_1Sa_positions_selftest_$(date +%Y%m%d_%H%M%S).log"
    ST="nad/results/08_1Sa_positions"
    log_info "NAD step 8 self-test: null"
    run_step8 "$ST/selftest_null"  "$ST/selftest_null"  --null
    log_info "NAD step 8 self-test: plant 30 %"
    run_step8 "$ST/selftest_plant" "$ST/selftest_plant" --plant 0.3
    log_info "self-test tables: $ST/selftest_{null,plant}/nad_signal.tsv"
    exit 0
fi

LOG="$NAD_DIR/logs/08_1Sa_positions_$(date +%Y%m%d_%H%M%S).log"
log_info "NAD step 8: NAD signal at CAGE single-A-clip positions"
run_step8 "nad/results/08_1Sa_positions" "nad/figures"
