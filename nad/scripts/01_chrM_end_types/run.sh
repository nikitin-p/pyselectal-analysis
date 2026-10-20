#!/usr/bin/env bash
# NAD branch step 1: the soft-clipped-A end type on the mitochondrial genome.
#
# Reads the Stage 8 typed CTSS table and the genome FASTA and annotation GTF
# named in the params file; no BAM is opened. Outputs go to
# nad/results/01_chrM_end_types/ and nad/figures/01_chrM_single_A_clip.{pdf,png}.
#
# Usage:
#   bash nad/scripts/01_chrM_end_types/run.sh [--params config/params.yaml] [--force]

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
LOG="$NAD_DIR/logs/01_chrM_end_types_$(date +%Y%m%d_%H%M%S).log"

log_info "NAD step 1: soft-clipped-A end type on chrM"
cd "$ROOT_DIR"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_chrM_end_types.R" \
    --ctss   "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
    --params "$PARAMS" \
    --outdir "nad/results/01_chrM_end_types" \
    --figdir "nad/figures" $FORCE_FLAG
