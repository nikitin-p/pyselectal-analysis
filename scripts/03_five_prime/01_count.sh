#!/usr/bin/env bash
# Stage 3, step 1: run pyselectal --count on every BAM in the sample manifest.
#
# Usage:
#   bash scripts/03_five_prime/01_count.sh [--samples FILE] [--outdir DIR] [--force]
#
# Defaults:
#   --samples  config/samples.tsv
#   --outdir   results/five_prime
#
# Columns required in samples.tsv: sample_id, bam
# (bam column may be a subsampled BAM path or point to results/bam_subsampled/)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
OUTDIR=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples) SAMPLES="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --force)   FORCE=1;      shift   ;;
        *) die "Unknown argument: $1" ;;
    esac
done

PYSELECTAL="$(yaml_get pyselectal "$PARAMS")"
CONDA_ENV="$(yaml_get conda_env_pyselectal "$PARAMS")"
COLLAPSE="$(yaml_get collapse_threshold "$PARAMS")"
MAPPED_PFX="$(yaml_get mapped_prefix "$PARAMS")"
[[ -z "$OUTDIR" ]] && OUTDIR="$(yaml_get results_five_prime "$PARAMS")"

# Resolve pyselectal path relative to project root
[[ "$PYSELECTAL" != /* ]] && PYSELECTAL="$ROOT_DIR/$PYSELECTAL"

require_cmd python
[[ -f "$PYSELECTAL" ]] || die "pyselectal not found: $PYSELECTAL"
[[ -f "$SAMPLES" ]]    || die "Samples file not found: $SAMPLES"

mkdir -p "$OUTDIR"

# Parse header to find column indices
header=$(grep -v '^#' "$SAMPLES" | head -1)
IFS=$'\t' read -ra cols <<< "$header"
sid_col=-1; bam_col=-1
for i in "${!cols[@]}"; do
    [[ "${cols[$i]}" == "sample_id" ]] && sid_col=$i
    [[ "${cols[$i]}" == "bam" ]]       && bam_col=$i
done
[[ $sid_col -lt 0 ]] && die "Column 'sample_id' not found in $SAMPLES"
[[ $bam_col -lt 0 ]] && die "Column 'bam' not found in $SAMPLES"

n_skipped=0; n_done=0; n_failed=0

while IFS=$'\t' read -ra row; do
    [[ "${row[0]}" =~ ^# ]] && continue
    [[ "${row[0]}" == "sample_id" ]] && continue
    sid="${row[$sid_col]}"
    bam="${row[$bam_col]}"

    # Resolve BAM path relative to project root if not absolute
    [[ "$bam" != /* ]] && bam="$ROOT_DIR/$bam"

    out="$OUTDIR/${sid}_counts.tsv"
    log="$ROOT_DIR/logs/03_five_prime_${sid}_$(date +%Y%m%d_%H%M%S).log"

    if [[ -f "$out" && $FORCE -eq 0 ]]; then
        log_info "[$sid] output exists, skipping (use --force to rerun)"
        (( n_skipped++ )) || true
        continue
    fi

    if [[ ! -f "$bam" ]]; then
        log_warn "[$sid] BAM not found: $bam — skipping"
        (( n_failed++ )) || true
        continue
    fi

    log_info "[$sid] counting 5' end types …"
    mkdir -p "$(dirname "$log")"
    if conda run -n "$CONDA_ENV" python "$PYSELECTAL" \
            -i "$bam" -c \
            --collapse-threshold "$COLLAPSE" \
            --mapped-prefix "$MAPPED_PFX" \
            -o "$out" \
            2>>"$log"; then
        log_info "[$sid] done → $out"
        (( n_done++ )) || true
    else
        log_error "[$sid] pyselectal failed — see $log"
        (( n_failed++ )) || true
    fi
done < <(grep -v '^#' "$SAMPLES")

log_info "Summary: done=$n_done  skipped=$n_skipped  failed=$n_failed"
[[ $n_failed -gt 0 ]] && exit 1 || exit 0
