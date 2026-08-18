#!/usr/bin/env bash
# Stage 2: subsample BAMs to equal depth.
# Usage: bash scripts/02_subsample/01_subsample.sh [--bam_dir DIR] [--out_dir DIR] [--force]
#
# 1. Count uniquely mapped reads (not alignments) in every *_uniq.bam.
# 2. Compute Q1, Q3, IQR across all libraries.
# 3. Exclude libraries whose depth < Q1 - 1.5×IQR (Tukey lower fence).
# 4. Set target = floor(subsample_factor * min(retained depths)).
# 5. Subsample each retained library to target using samtools view -s SEED.FRAC.
# 6. Write the depth table and the summary statistics, then draw the depth plot.
#
# Skips samples whose output already exists unless --force is given.
# Compatible with bash 3.2+ (macOS default).

set -euo pipefail
export LC_NUMERIC=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_ROOT/scripts/utils/log.sh"
source "$PROJECT_ROOT/scripts/utils/parse_yaml.sh"

FORCE=0
BAM_IN_DIR=""
BAM_OUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam_dir) BAM_IN_DIR="$2";  shift 2 ;;
        --params)  PARAMS_ARG="$2";  shift 2 ;;
        --out_dir) BAM_OUT_DIR="$2"; shift 2 ;;
        --force)   FORCE=1;          shift   ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# ── config ────────────────────────────────────────────────────────────────────
PARAMS="${PARAMS_ARG:-$PROJECT_ROOT/config/params.yaml}"
[[ -z "$BAM_IN_DIR"  ]] && BAM_IN_DIR="$(yaml_get results_bam     "$PARAMS")"
[[ -z "$BAM_OUT_DIR" ]] && BAM_OUT_DIR="$(yaml_get results_bam_sub "$PARAMS")"
[[ "$BAM_IN_DIR"  != /* ]] && BAM_IN_DIR="$PROJECT_ROOT/$BAM_IN_DIR"
[[ "$BAM_OUT_DIR" != /* ]] && BAM_OUT_DIR="$PROJECT_ROOT/$BAM_OUT_DIR"
LOG_DIR="$PROJECT_ROOT/logs"

RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$PROJECT_ROOT")"
QC_DIR="$(yaml_get results_qc "$PARAMS")"
[[ -z "$QC_DIR" ]] && QC_DIR="$RESULTS_ROOT/qc"
[[ "$QC_DIR" != /* ]] && QC_DIR="$PROJECT_ROOT/$QC_DIR"
FIG_DIR="$(yaml_get results_figures "$PARAMS")"
[[ -z "$FIG_DIR" ]] && FIG_DIR="$RESULTS_ROOT/figures"
[[ "$FIG_DIR" != /* ]] && FIG_DIR="$PROJECT_ROOT/$FIG_DIR"

DEPTH_TSV="$QC_DIR/library_depth.tsv"
STATS_TSV="$QC_DIR/depth_stats.tsv"

CONDA_TOOLS="$(yaml_get conda_env_tools "$PARAMS")"
CONDA_R="$(yaml_get conda_env_r "$PARAMS")"
SAMTOOLS=$(conda run -n "$CONDA_TOOLS" which samtools)

SEED="$(yaml_get subsample_seed "$PARAMS")"
FACTOR="$(yaml_get subsample_factor "$PARAMS")"
THREADS=4
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$BAM_OUT_DIR" "$LOG_DIR" "$QC_DIR" "$FIG_DIR"

# ── collect BAMs and count reads (using parallel arrays for bash 3.2 compat) ─
SAMPLES=()
COUNTS=()
BAMS=()

for BAM in "$BAM_IN_DIR"/*_uniq.bam; do
    [[ -f "$BAM" ]] || continue
    BAMS+=("$BAM")
done

[[ ${#BAMS[@]} -eq 0 ]] && die "No *_uniq.bam files found in $BAM_IN_DIR"

log_info "Found ${#BAMS[@]} BAMs — counting reads..."

for BAM in "${BAMS[@]}"; do
    SAMPLE=$(basename "$BAM" _uniq.bam)
    N=$($SAMTOOLS view -c -F 256 "$BAM")
    SAMPLES+=("$SAMPLE")
    COUNTS+=("$N")
    log_info "  $SAMPLE: $N"
done

# ── compute Q1, Q3, IQR and Tukey lower fence ──────────────────────────────────
read -r Q1 Q3 IQR DEPTH_THRESHOLD < <(printf '%s\n' "${COUNTS[@]}" | sort -n | awk '
    {a[NR]=$1}
    END{
        n=NR
        q1_pos = (n+1)*0.25
        q3_pos = (n+1)*0.75
        q1_lo = int(q1_pos); q1_hi = q1_lo + 1
        q3_lo = int(q3_pos); q3_hi = q3_lo + 1
        q1_frac = q1_pos - q1_lo
        q3_frac = q3_pos - q3_lo
        if (q1_hi > n) q1_hi = n
        if (q3_hi > n) q3_hi = n
        Q1 = a[q1_lo] + q1_frac * (a[q1_hi] - a[q1_lo])
        Q3 = a[q3_lo] + q3_frac * (a[q3_hi] - a[q3_lo])
        IQR = Q3 - Q1
        threshold = Q1 - 1.5 * IQR
        if (threshold < 0) threshold = 0
        printf "%d %d %d %d\n", Q1, Q3, IQR, threshold
    }')
log_info "Q1=$Q1  Q3=$Q3  IQR=$IQR  Tukey lower fence (exclusion threshold)=$DEPTH_THRESHOLD"

# ── filter low-coverage libraries ─────────────────────────────────────────────
RETAINED_SAMPLES=()
RETAINED_COUNTS=()
EXCLUDED=0

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    N="${COUNTS[$i]}"
    if (( N < DEPTH_THRESHOLD )); then
        log_warn "EXCLUDED $SAMPLE: $N reads < threshold $DEPTH_THRESHOLD (Q1 - 1.5×IQR)"
        EXCLUDED=$((EXCLUDED + 1))
    else
        RETAINED_SAMPLES+=("$SAMPLE")
        RETAINED_COUNTS+=("$N")
        log_info "RETAINED $SAMPLE: $N reads"
    fi
done

[[ ${#RETAINED_SAMPLES[@]} -eq 0 ]] && die "All libraries excluded — check min_library_depth_fraction in params.yaml"

# ── compute subsampling target from retained libraries ────────────────────────
MIN_RETAINED=999999999999
for N in "${RETAINED_COUNTS[@]}"; do
    (( N < MIN_RETAINED )) && MIN_RETAINED=$N
done

TARGET=$(echo "$MIN_RETAINED $FACTOR" | awk '{printf "%d", $1 * $2}')
log_info "Retained ${#RETAINED_SAMPLES[@]} / ${#SAMPLES[@]} libraries — min=$MIN_RETAINED  factor=$FACTOR  target=$TARGET"

# ── depth table and summary statistics ────────────────────────────────────────
# These are written on every run, including runs in which every subsampled BAM
# already exists, because the counting loop above is unconditional and the two
# tables are the only durable record of the depths the exclusion rule saw.
{
    printf 'sample_id\tdepth\tretained\n'
    for i in "${!SAMPLES[@]}"; do
        N="${COUNTS[$i]}"
        if (( N < DEPTH_THRESHOLD )); then KEEP="no"; else KEEP="yes"; fi
        printf '%s\t%s\t%s\n' "${SAMPLES[$i]}" "$N" "$KEEP"
    done
} > "$DEPTH_TSV"

{
    printf 'metric\tvalue\n'
    printf 'n_libraries\t%s\n'       "${#SAMPLES[@]}"
    printf 'n_retained\t%s\n'        "${#RETAINED_SAMPLES[@]}"
    printf 'n_excluded\t%s\n'        "$EXCLUDED"
    printf 'q1\t%s\n'                "$Q1"
    printf 'q3\t%s\n'                "$Q3"
    printf 'iqr\t%s\n'               "$IQR"
    printf 'tukey_lower_fence\t%s\n' "$DEPTH_THRESHOLD"
    printf 'min_retained\t%s\n'      "$MIN_RETAINED"
    printf 'subsample_factor\t%s\n'  "$FACTOR"
    printf 'subsample_target\t%s\n'  "$TARGET"
} > "$STATS_TSV"

log_info "Depth table → $DEPTH_TSV"
log_info "Depth statistics → $STATS_TSV"

# ── subsample ─────────────────────────────────────────────────────────────────
for i in "${!RETAINED_SAMPLES[@]}"; do
    SAMPLE="${RETAINED_SAMPLES[$i]}"
    N="${RETAINED_COUNTS[$i]}"
    BAM="$BAM_IN_DIR/${SAMPLE}_uniq.bam"
    OUT="$BAM_OUT_DIR/${SAMPLE}_sub.bam"
    LOG="$LOG_DIR/02_subsample_${SAMPLE}_$(date '+%Y%m%d_%H%M%S').log"

    if [[ -f "$OUT" && $FORCE -eq 0 ]]; then
        log_info "$SAMPLE — output exists, skipping (use --force to overwrite)"
        continue
    fi

    FRAC=$(echo "$N $TARGET" | awk '{f=$2/$1; if(f>1)f=1; printf "%.6f", f}')
    FRAC_PART="${FRAC#0}"
    S_FLAG="${SEED}${FRAC_PART}"

    log_info "$SAMPLE — N=$N  frac=$FRAC  -s $S_FLAG"

    $SAMTOOLS view -b -s "$S_FLAG" -@ "$THREADS" "$BAM" \
        | $SAMTOOLS sort -@ "$THREADS" -o "$OUT" \
        2> "$LOG"

    $SAMTOOLS index "$OUT"
    ACTUAL=$($SAMTOOLS view -c "$OUT")
    log_info "$SAMPLE — done: $ACTUAL reads → $OUT"
done

# ── depth plot ────────────────────────────────────────────────────────────────
# The plot is a report of the exclusion decision rather than an input to it, so a
# failure here leaves the subsampled BAMs and both tables intact and is reported
# as a warning rather than allowed to fail the stage.
PLOT_R="$SCRIPT_DIR/02_plot_depth.R"
PLOT_LOG="$LOG_DIR/02_subsample_depth_plot_$(date '+%Y%m%d_%H%M%S').log"
if RSCRIPT=$(conda run -n "$CONDA_R" which Rscript 2>/dev/null) && [[ -n "$RSCRIPT" ]]; then
    if "$RSCRIPT" "$PLOT_R" \
            --depth  "$DEPTH_TSV" \
            --stats  "$STATS_TSV" \
            --outdir "$FIG_DIR" \
            > "$PLOT_LOG" 2>&1; then
        log_info "Depth plot → $FIG_DIR/02_library_depth.pdf"
    else
        log_warn "Depth plot failed — see $PLOT_LOG"
    fi
else
    log_warn "Rscript not found in conda env '$CONDA_R' — skipping the depth plot"
fi

log_info "Stage 2 complete."
