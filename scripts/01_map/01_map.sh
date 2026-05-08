#!/usr/bin/env bash
# Stage 1: map FASTQ files with STAR and filter for uniquely mapped reads (MAPQ = 255).
#
# Usage:
#   bash scripts/01_map/01_map.sh [--samples FILE] [--outdir DIR] [--force]
#
# Defaults:
#   --samples  config/samples.tsv
#   --outdir   results/bam   (from params.yaml results_bam)
#
# Reads from samples.tsv: sample_id, dataset, fastq.
# Skips rows with an empty fastq column or no fastq path.
# STAR index is selected by dataset:
#   yeast_cage              -> star_index_saccer3
#   luhmes_cage / luhmes_netcage / netcage_* -> star_index_hg38
#
# Outputs: {outdir}/{sample_id}_uniq.bam  (sorted, indexed, MAPQ-filtered)
# Skips samples whose output already exists unless --force is given.

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

# ── read config ───────────────────────────────────────────────────────────────
THREADS="$(yaml_get threads "$PARAMS")"
MAPQ="$(yaml_get mapq_threshold "$PARAMS")"
STAR_EXTRA_RAW="$(yaml_get star_extra_args "$PARAMS")"
STAR_IDX_SACCER3="$(yaml_get star_index_saccer3 "$PARAMS")"
STAR_IDX_HG38="$(yaml_get star_index_hg38 "$PARAMS")"
CONDA_STAR="$(yaml_get conda_env_star "$PARAMS")"
CONDA_TOOLS="$(yaml_get conda_env_tools "$PARAMS")"
[[ -z "$OUTDIR" ]] && OUTDIR="$ROOT_DIR/$(yaml_get results_bam "$PARAMS")"

# Split extra args string into array for safe quoting
read -ra STAR_EXTRA <<< "$STAR_EXTRA_RAW"

STAR=$(conda run -n "$CONDA_STAR" which STAR)
SAMTOOLS=$(conda run -n "$CONDA_TOOLS" which samtools)
[[ -f "$SAMPLES" ]] || die "Samples file not found: $SAMPLES"

mkdir -p "$OUTDIR" "$ROOT_DIR/logs"

# ── parse header ──────────────────────────────────────────────────────────────
header=$(grep -v '^#' "$SAMPLES" | head -1)
IFS=$'\t' read -ra cols <<< "$header"
sid_col=-1; dataset_col=-1; fastq_col=-1
for i in "${!cols[@]}"; do
    [[ "${cols[$i]}" == "sample_id" ]] && sid_col=$i
    [[ "${cols[$i]}" == "dataset"   ]] && dataset_col=$i
    [[ "${cols[$i]}" == "fastq"     ]] && fastq_col=$i
done
[[ $sid_col -lt 0     ]] && die "Column 'sample_id' not found in $SAMPLES"
[[ $dataset_col -lt 0 ]] && die "Column 'dataset' not found in $SAMPLES"
[[ $fastq_col -lt 0   ]] && die "Column 'fastq' not found in $SAMPLES"

n_skipped=0; n_done=0; n_failed=0

# ── process samples ───────────────────────────────────────────────────────────
while IFS=$'\t' read -ra row; do
    [[ "${row[0]}" =~ ^# ]] && continue
    [[ "${row[0]}" == "sample_id" ]] && continue

    sid="${row[$sid_col]}"
    dataset="${row[$dataset_col]}"
    fastq="${row[$fastq_col]:-}"

    [[ -z "$fastq" ]] && continue   # skip rows with no FASTQ (e.g. BAM-only entries)

    # Resolve FASTQ path relative to project root
    [[ "$fastq" != /* ]] && fastq="$ROOT_DIR/$fastq"

    # Select STAR index by dataset
    case "$dataset" in
        yeast_cage)
            STAR_IDX="$STAR_IDX_SACCER3"
            [[ -z "$STAR_IDX" ]] && die "[$sid] star_index_saccer3 not set in params.yaml"
            ;;
        luhmes_cage|luhmes_netcage|netcage_*)
            STAR_IDX="$STAR_IDX_HG38"
            [[ -z "$STAR_IDX" ]] && die "[$sid] star_index_hg38 not set in params.yaml"
            ;;
        *)
            die "[$sid] Unknown dataset '$dataset' — add a case to select a STAR index"
            ;;
    esac

    OUT_BAM="$OUTDIR/${sid}_uniq.bam"
    TMP_DIR="$OUTDIR/${sid}_star_tmp"
    LOG="$ROOT_DIR/logs/01_map_${sid}_$(date '+%Y%m%d_%H%M%S').log"

    if [[ -f "$OUT_BAM" && $FORCE -eq 0 ]]; then
        log_info "[$sid] output exists, skipping (use --force to rerun)"
        (( n_skipped++ )) || true
        continue
    fi

    if [[ ! -f "$fastq" ]]; then
        log_warn "[$sid] FASTQ not found: $fastq — skipping"
        (( n_failed++ )) || true
        continue
    fi

    log_info "[$sid] mapping with STAR → $TMP_DIR"
    log_info "[$sid]   index : $STAR_IDX"
    log_info "[$sid]   fastq : $fastq"
    mkdir -p "$TMP_DIR"

    $STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_IDX" \
        --readFilesIn "$fastq" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 \
        --outFilterMultimapNmax 1 \
        --limitBAMsortRAM 3000000000 \
        --outFileNamePrefix "$TMP_DIR/" \
        "${STAR_EXTRA[@]}" \
        > "$LOG" 2>&1

    RAW_BAM="$TMP_DIR/Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$RAW_BAM" ]]; then
        log_error "[$sid] STAR output BAM not found — see $LOG"
        (( n_failed++ )) || true
        continue
    fi

    log_info "[$sid] filtering MAPQ ≥ $MAPQ …"
    $SAMTOOLS view -b -q "$MAPQ" -@ "$THREADS" "$RAW_BAM" -o "$OUT_BAM" 2>> "$LOG"
    $SAMTOOLS index "$OUT_BAM"

    MAPPED=$($SAMTOOLS view -c "$OUT_BAM")
    log_info "[$sid] done — $MAPPED uniquely mapped reads → $OUT_BAM"

    rm -rf "$TMP_DIR"
    (( n_done++ )) || true

done < <(grep -v '^#' "$SAMPLES")

log_info "Summary: done=$n_done  skipped=$n_skipped  failed=$n_failed"
[[ $n_failed -gt 0 ]] && exit 1 || exit 0
