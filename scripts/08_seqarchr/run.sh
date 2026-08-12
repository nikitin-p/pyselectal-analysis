#!/usr/bin/env bash
# Stage 08: de novo sequence architectures (seqArchR) by 5' end type
#
# Asks whether the genomic context around TSSs partitions into distinct
# sequence architectures, and whether that partition is associated with 5' end
# type. seqArchR never sees the end-type labels, so any association it produces
# is an independent result rather than a circular one.
#
# Usage:
#   bash scripts/08_seqarchr/run.sh [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
FORCE_FLAG=""
[[ "${1:-}" == "--force" ]] && FORCE_FLAG="--force"

OUTDIR="$ROOT_DIR/results/seqarchr"
FIGDIR="$ROOT_DIR/results/figures"

UPSTREAM="$(yaml_get seqarchr_upstream   "$PARAMS")"
DOWNSTREAM="$(yaml_get seqarchr_downstream "$PARAMS")"
MIN_TOTAL="$(yaml_get seqarchr_min_total  "$PARAMS")"
TOTAL_ITR="$(yaml_get seqarchr_total_itr  "$PARAMS")"
CHUNK_SIZE="$(yaml_get seqarchr_chunk_size "$PARAMS")"
CV_FOLDS="$(yaml_get seqarchr_cv_folds    "$PARAMS")"
N_CORES="$(yaml_get seqarchr_n_cores      "$PARAMS")"
N_ARCH="$(yaml_get seqarchr_n_arch        "$PARAMS")"
SEED="$(yaml_get seqarchr_seed            "$PARAMS")"
PYENV="$(yaml_get conda_env_seqarchr      "$PARAMS")"

# seqArchR factorises via scikit-learn through reticulate, so the Python
# interpreter must be pinned to the env that carries scikit-learn.
PYBIN="$(conda run -n "$PYENV" which python)"
export RETICULATE_PYTHON="$PYBIN"
log_info "Python backend: $RETICULATE_PYTHON"

# seqArchR is installed in the system R library, not in conda_env_r.
RSCRIPT="$(which Rscript)"

log_info "Stage 08 step 1: extracting promoter windows"
$RSCRIPT "$SCRIPT_DIR/01_prepare_seqs.R" \
    --tss        "$ROOT_DIR/results/tss/tss_clustered.tsv" \
    --genome     "$ROOT_DIR/data/sacCer3.fa" \
    --outdir     "$OUTDIR" \
    --upstream   "$UPSTREAM" \
    --downstream "$DOWNSTREAM" \
    --min_total  "$MIN_TOTAL" \
    $FORCE_FLAG

log_info "Stage 08 step 2: running seqArchR"
$RSCRIPT "$SCRIPT_DIR/02_run_seqarchr.R" \
    --fasta      "$OUTDIR/promoter_windows.fa" \
    --outdir     "$OUTDIR" \
    --upstream   "$UPSTREAM" \
    --downstream "$DOWNSTREAM" \
    --total_itr  "$TOTAL_ITR" \
    --chunk_size "$CHUNK_SIZE" \
    --cv_folds   "$CV_FOLDS" \
    --n_cores    "$N_CORES" \
    --n_arch     "$N_ARCH" \
    --seed       "$SEED" \
    $FORCE_FLAG

log_info "Stage 08 step 3: architecture vs end-type analysis"
$RSCRIPT "$SCRIPT_DIR/03_analyze_clusters.R" \
    --outdir     "$OUTDIR" \
    --figdir     "$FIGDIR" \
    --upstream   "$UPSTREAM" \
    --downstream "$DOWNSTREAM" \
    $FORCE_FLAG

log_info "Stage 08 complete."
log_info "  tables  → results/seqarchr/"
log_info "  figures → results/figures/16_*.pdf"
