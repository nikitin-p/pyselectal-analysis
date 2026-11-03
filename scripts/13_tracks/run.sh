#!/usr/bin/env bash
# Stage 21: stranded CTSS bigWigs from the subsampled libraries.
#
# One black track holds the 5' ends pooled over every subsampled library, and
# one coloured track per 5' end type holds every read of that type, so mixed
# and dominated positions alike are read off the hub by overlaying the tracks.
# track_manifest.tsv lists each file with its colour and totals.
#
# This stage reads Stage 17's fine-typed CTSS table only, so no BAM is opened.
# Run scripts/10_dominant_endtype/run.sh first.
#
# Usage:
#   bash scripts/13_tracks/run.sh [--params ...] [--samples ...] [--force]

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

resolve_path() {
    local v="$1"
    [[ -z "$v" ]] && { printf ''; return; }
    [[ "$v" != /* ]] && v="$ROOT_DIR/$v"
    printf '%s' "$v"
}
GENOME_FASTA="$(resolve_path "$(yaml_get genome_fasta "$PARAMS")")"
DOM_DIR="$ROOT_DIR/$(yaml_get results_dominance "$PARAMS")"
TRACK_DIR="$ROOT_DIR/$(yaml_get results_tracks_sub "$PARAMS")"

LOG="$ROOT_DIR/logs/21_tracks_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 21: subsampled 5' end bigWigs"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_endtype_bigwigs.R" \
    --ctss      "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --fai       "$GENOME_FASTA.fai" \
    --outdir    "$TRACK_DIR" $FORCE_FLAG

# The hub reads all three manifests, so it also picks up the Stage 10 and
# Stage 22 tracks when they exist and leaves them out otherwise.
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_make_hub.R" \
    --endtype  "$TRACK_DIR/track_manifest.tsv" \
    --cager    "$ROOT_DIR/$(yaml_get results_cager "$PARAMS")/tracks/track_manifest.tsv" \
    --longclip "$ROOT_DIR/$(yaml_get results_longclip "$PARAMS")/tracks/track_manifest.tsv" \
    --out      "$TRACK_DIR/hub.txt"

log_info "Stage 21 complete."
log_info "  hub      → $TRACK_DIR/hub.txt"
log_info "  bigWigs  → $TRACK_DIR"
log_info "  manifest → $TRACK_DIR/track_manifest.tsv"
