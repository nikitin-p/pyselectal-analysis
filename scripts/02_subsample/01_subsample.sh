#!/usr/bin/env bash
# Stage 2: subsample BAMs to equal depth.
# Usage: bash scripts/02_subsample/01_subsample.sh [--force]
#
# 1. Count reads in every *_uniq.bam.
# 2. Compute the median depth across all libraries.
# 3. Exclude libraries whose depth < median * min_library_depth_fraction (from params.yaml).
# 4. Set target = floor(subsample_factor * min(retained depths)).
# 5. Subsample each retained library to target using samtools view -s SEED.FRAC.
#
# Skips samples whose output already exists unless --force is given.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_ROOT/scripts/utils/log.sh"
source "$PROJECT_ROOT/scripts/utils/parse_yaml.sh"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

# ── config ────────────────────────────────────────────────────────────────────
PARAMS="$PROJECT_ROOT/config/params.yaml"
BAM_IN_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/results/bam"
BAM_OUT_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/results/bam_subsampled"
LOG_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/logs"

SEED="$(yaml_get subsample_seed "$PARAMS")"
FACTOR="$(yaml_get subsample_factor "$PARAMS")"
MIN_FRAC="$(yaml_get min_library_depth_fraction "$PARAMS")"
THREADS=4
# ─────────────────────────────────────────────────────────────────────────────

require_cmd samtools
mkdir -p "$BAM_OUT_DIR" "$LOG_DIR"

# ── collect BAMs and count reads ──────────────────────────────────────────────
mapfile -t BAMS < <(ls "$BAM_IN_DIR"/*_uniq.bam 2>/dev/null)
[[ ${#BAMS[@]} -eq 0 ]] && die "No *_uniq.bam files found in $BAM_IN_DIR"

log_info "Found ${#BAMS[@]} BAMs — counting reads..."

declare -A COUNTS
for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" _uniq.bam)
    N=$(samtools view -c "$BAM")
    COUNTS[$SAMPLE]=$N
    log_info "  $SAMPLE: $N"
done

# ── compute median depth ───────────────────────────────────────────────────────
ALL_COUNTS=( "${COUNTS[@]}" )
MEDIAN=$(printf '%s\n' "${ALL_COUNTS[@]}" | sort -n | awk '{a[NR]=$1}
    END{
        n=NR
        if(n%2==1) print a[(n+1)/2]
        else print int((a[n/2]+a[n/2+1])/2)
    }')
DEPTH_THRESHOLD=$(echo "$MEDIAN $MIN_FRAC" | awk '{printf "%d", $1 * $2}')
log_info "Median depth=$MEDIAN  min_fraction=$MIN_FRAC  exclusion threshold=$DEPTH_THRESHOLD"

# ── filter low-coverage libraries ─────────────────────────────────────────────
declare -A RETAINED
for SAMPLE in "${!COUNTS[@]}"; do
    N=${COUNTS[$SAMPLE]}
    if (( N < DEPTH_THRESHOLD )); then
        log_warn "EXCLUDED $SAMPLE: $N reads < threshold $DEPTH_THRESHOLD (< $MIN_FRAC × median $MEDIAN)"
    else
        RETAINED[$SAMPLE]=$N
        log_info "RETAINED $SAMPLE: $N reads"
    fi
done

[[ ${#RETAINED[@]} -eq 0 ]] && die "All libraries excluded — check min_library_depth_fraction in params.yaml"

# ── compute subsampling target from retained libraries ────────────────────────
MIN_RETAINED=999999999999
for N in "${RETAINED[@]}"; do
    (( N < MIN_RETAINED )) && MIN_RETAINED=$N
done

TARGET=$(echo "$MIN_RETAINED $FACTOR" | awk '{printf "%d", $1 * $2}')
log_info "Retained ${#RETAINED[@]} / ${#COUNTS[@]} libraries — min=$MIN_RETAINED  factor=$FACTOR  target=$TARGET"

# ── subsample ─────────────────────────────────────────────────────────────────
for SAMPLE in "${!RETAINED[@]}"; do
    BAM="$BAM_IN_DIR/${SAMPLE}_uniq.bam"
    OUT="$BAM_OUT_DIR/${SAMPLE}_sub.bam"
    LOG="$LOG_DIR/02_subsample_${SAMPLE}_$(date '+%Y%m%d_%H%M%S').log"

    if [[ -f "$OUT" && $FORCE -eq 0 ]]; then
        log_info "$SAMPLE — output exists, skipping (use --force to overwrite)"
        continue
    fi

    N=${RETAINED[$SAMPLE]}
    FRAC=$(echo "$N $TARGET" | awk '{f=$2/$1; if(f>1)f=1; printf "%.6f", f}')
    FRAC_PART="${FRAC#0}"
    S_FLAG="${SEED}${FRAC_PART}"

    log_info "$SAMPLE — N=$N  frac=$FRAC  -s $S_FLAG"

    samtools view -b -s "$S_FLAG" -@ "$THREADS" "$BAM" \
        | samtools sort -@ "$THREADS" -o "$OUT" \
        2> "$LOG"

    samtools index "$OUT"
    ACTUAL=$(samtools view -c "$OUT")
    log_info "$SAMPLE — done: $ACTUAL reads → $OUT"
done

log_info "Stage 2 complete."
