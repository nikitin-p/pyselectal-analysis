#!/usr/bin/env bash
# Stage 2: subsample BAMs to equal depth.
# Usage: bash scripts/02_subsample/01_subsample.sh [--force]
#
# Reads counts from BAM files, computes target = floor(0.9 * min_count),
# then runs samtools view -s SEED.FRAC for each sample.
# Skips samples whose output already exists unless --force is given.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_ROOT/scripts/utils/log.sh"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

# ── config ────────────────────────────────────────────────────────────────────
BAM_IN_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/results/bam"
BAM_OUT_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/results/bam_subsampled"
LOG_DIR="/home/pnikitin/cage/analysis_pyselectal_for_paper/logs"
SEED=42
FACTOR=0.9
THREADS=4
# ─────────────────────────────────────────────────────────────────────────────

require_cmd samtools
mkdir -p "$BAM_OUT_DIR" "$LOG_DIR"

# Collect all *_uniq.bam files
mapfile -t BAMS < <(ls "$BAM_IN_DIR"/*_uniq.bam 2>/dev/null)
[[ ${#BAMS[@]} -eq 0 ]] && die "No *_uniq.bam files found in $BAM_IN_DIR"

log_info "Found ${#BAMS[@]} BAMs — counting reads..."

declare -A COUNTS
MIN_COUNT=999999999999

for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" _uniq.bam)
    N=$(samtools view -c "$BAM")
    COUNTS[$SAMPLE]=$N
    log_info "  $SAMPLE: $N"
    (( N < MIN_COUNT )) && MIN_COUNT=$N
done

TARGET=$(echo "$MIN_COUNT $FACTOR" | awk '{printf "%d", $1 * $2}')
log_info "min=$MIN_COUNT  factor=$FACTOR  target=$TARGET"

# ── subsample ─────────────────────────────────────────────────────────────────
for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" _uniq.bam)
    OUT="$BAM_OUT_DIR/${SAMPLE}_sub.bam"
    LOG="$LOG_DIR/02_subsample_${SAMPLE}_$(date '+%Y%m%d_%H%M%S').log"

    if [[ -f "$OUT" && $FORCE -eq 0 ]]; then
        log_info "$SAMPLE — output exists, skipping (use --force to overwrite)"
        continue
    fi

    N=${COUNTS[$SAMPLE]}
    # Compute fraction to 6 decimal places; cap at 1.0
    FRAC=$(echo "$N $TARGET" | awk '{f=$2/$1; if(f>1)f=1; printf "%.6f", f}')
    # Strip leading zero to get the decimal part for -s INT.FRAC format
    FRAC_PART="${FRAC#0}"   # e.g. .258384
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
