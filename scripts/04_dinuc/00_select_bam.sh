#!/usr/bin/env bash
# Stage 6 pre-step: filter subsampled BAMs by 5' end type using pyselectal --select.
# Produces one BAM per sample per type, and a per-type samples TSV for dinuc analysis.
#
# Usage:
#   bash scripts/04_dinuc/00_select_bam.sh [--bam_dir DIR] [--out_dir DIR] \
#       [--samples FILE] [--params FILE] [--types "1Sg 2Sg"] [--force]
#
# Outputs:
#   <out_dir>/{sample_id}_1Sg.bam  (and .bai)
#   <out_dir>/{sample_id}_2Sg.bam  (and .bai)
#   config/samples_1Sg.tsv
#   config/samples_2Sg.tsv

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
BAM_DIR=""
OUT_DIR=""
TYPES="1Sg 2Sg"
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam_dir)  BAM_DIR="$2";  shift 2 ;;
        --out_dir)  OUT_DIR="$2";  shift 2 ;;
        --samples)  SAMPLES="$2";  shift 2 ;;
        --params)   PARAMS="$2";   shift 2 ;;
        --types)    TYPES="$2";    shift 2 ;;
        --force)    FORCE=1;       shift   ;;
        *) die "Unknown argument: $1" ;;
    esac
done

PYSELECTAL="$(yaml_get pyselectal "$PARAMS")"
CONDA_ENV="$(yaml_get conda_env_pyselectal "$PARAMS")"
[[ "$PYSELECTAL" != /* ]] && PYSELECTAL="$ROOT_DIR/$PYSELECTAL"
[[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$BAM_DIR")/bam_selected"

require_cmd samtools
[[ -f "$PYSELECTAL" ]] || die "pyselectal not found: $PYSELECTAL"
[[ -f "$SAMPLES" ]]    || die "Samples file not found: $SAMPLES"
[[ -n "$BAM_DIR" ]]    || die "--bam_dir is required"

mkdir -p "$OUT_DIR"

# pyselectal --select spec for each type label
select_spec() {
    case "$1" in
        1Sg) echo "1Sg"  ;;
        2Sg) echo "2.Sg" ;;
        *)   echo "$1"   ;;
    esac
}

# Parse header
header=$(grep -v '^#' "$SAMPLES" | head -1)
IFS=$'\t' read -ra cols <<< "$header"
sid_col=-1; bam_col=-1
for i in "${!cols[@]}"; do
    [[ "${cols[$i]}" == "sample_id" ]] && sid_col=$i
    [[ "${cols[$i]}" == "bam" ]]       && bam_col=$i
done
[[ $sid_col -lt 0 ]] && die "Column 'sample_id' not found"
[[ $bam_col -lt 0 ]] && die "Column 'bam' not found"

# Build per-type samples TSVs (header line first)
for TYPE in $TYPES; do
    TSV_OUT="$ROOT_DIR/config/samples_${TYPE}.tsv"
    echo "$header" > "$TSV_OUT"
done

# Process each sample
while IFS=$'\t' read -ra row; do
    [[ "${row[0]}" =~ ^# ]] && continue
    [[ "${row[0]}" == "sample_id" ]] && continue
    sid="${row[$sid_col]}"
    bam="${row[$bam_col]}"
    [[ "$bam" != /* ]] && bam="$ROOT_DIR/$bam"

    for TYPE in $TYPES; do
        SPEC="$(select_spec "$TYPE")"
        OUT_BAM="$OUT_DIR/${sid}_${TYPE}.bam"

        if [[ -f "$OUT_BAM" && $FORCE -eq 0 ]]; then
            log_info "[$sid/$TYPE] exists, skipping"
        else
            log_info "[$sid/$TYPE] selecting spec=$SPEC from $bam"
            conda run -n "$CONDA_ENV" python "$PYSELECTAL" \
                -i "$bam" --select "$SPEC" -o "$OUT_BAM"
            samtools index "$OUT_BAM"
            log_info "[$sid/$TYPE] done → $OUT_BAM"
        fi

        # Append row to per-type TSV with updated bam column
        row[$bam_col]="$OUT_BAM"
        printf '%s' "${row[0]}"; for c in "${row[@]:1}"; do printf '\t%s' "$c"; done; printf '\n'
    done >> /dev/null  # rows written inside loop below

done < <(grep -v '^#' "$SAMPLES")

# Re-build per-type TSVs cleanly (done separately to avoid subshell issues)
for TYPE in $TYPES; do
    SPEC="$(select_spec "$TYPE")"
    TSV_OUT="$ROOT_DIR/config/samples_${TYPE}.tsv"
    # header
    echo "$header" > "$TSV_OUT"
    # rows: replace bam column with selected BAM path
    while IFS=$'\t' read -ra row; do
        [[ "${row[0]}" =~ ^# ]] && continue
        [[ "${row[0]}" == "sample_id" ]] && continue
        sid="${row[$sid_col]}"
        row[$bam_col]="$OUT_DIR/${sid}_${TYPE}.bam"
        printf '%s' "${row[0]}"
        for c in "${row[@]:1}"; do printf '\t%s' "$c"; done
        printf '\n'
    done < <(grep -v '^#' "$SAMPLES") >> "$TSV_OUT"
    log_info "Written: $TSV_OUT"
done

log_info "Selection complete. Run dinuc with --samples config/samples_1Sg.tsv or samples_2Sg.tsv"
