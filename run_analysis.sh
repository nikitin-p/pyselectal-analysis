#!/usr/bin/env bash
# Top-level runner: stages 1–9 (mapping → dinucleotides → TSS matrix).
#
# Usage:
#   bash run_analysis.sh [--force]
#
# Prerequisites:
#   - FASTQ paths filled in config/samples.tsv
#   - STAR indexes built (star_index_saccer3 / star_index_hg38 in params.yaml)

set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
[[ "${1:-}" == "--force" ]] && FORCE_FLAG="--force"

CONDA_R="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_R" which Rscript)

BAM_SUB="$ROOT_DIR/results/bam_subsampled"
BAM_SEL="$ROOT_DIR/results/bam_selected"

# ── Stage 1: STAR mapping ─────────────────────────────────────────────────────
log_info "Stage 1: STAR mapping"
bash "$ROOT_DIR/scripts/01_map/01_map.sh" \
    --samples "$SAMPLES" \
    $FORCE_FLAG

# ── Stage 2: subsampling ──────────────────────────────────────────────────────
log_info "Stage 2: subsampling to equal depth"
bash "$ROOT_DIR/scripts/02_subsample/01_subsample.sh" $FORCE_FLAG

# ── Stage 3: pyselectal --count ───────────────────────────────────────────────
log_info "Stage 3: 5' end type counts"
bash "$ROOT_DIR/scripts/03_five_prime/01_count.sh" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/results/five_prime" \
    $FORCE_FLAG

# ── Stage 4: per-sample histograms + heatmap + chi-square ────────────────────
log_info "Stage 4: replicate heatmap and chi-square"
$RSCRIPT "$ROOT_DIR/scripts/03_five_prime/02_plot.R" \
    --counts_dir "$ROOT_DIR/results/five_prime" \
    --samples    "$SAMPLES" \
    --outdir     "$ROOT_DIR/results/figures" \
    --params     "$PARAMS"

# ── Stage 5: soft-clip nucleotide composition ─────────────────────────────────
log_info "Stage 5: soft-clip nucleotide composition"
$RSCRIPT "$ROOT_DIR/scripts/03_five_prime/03_softclip_composition.R" \
    --counts_dir "$ROOT_DIR/results/five_prime" \
    --samples    "$SAMPLES" \
    --outdir     "$ROOT_DIR/results/figures" \
    --params     "$PARAMS"

# ── Stage 6a: select per-type BAMs ───────────────────────────────────────────
log_info "Stage 6a: select 1Sg and 2Sg BAMs"
bash "$ROOT_DIR/scripts/04_dinuc/00_select_bam.sh" \
    --bam_dir "$BAM_SUB" \
    $FORCE_FLAG

# ── Stage 6b: dinucleotide proportions (1Sg) ─────────────────────────────────
log_info "Stage 6b: dinucleotide proportions — 1Sg"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SEL" \
    --samples "$ROOT_DIR/config/samples_1Sg.tsv" \
    --outdir  "$ROOT_DIR/results/dinuc_1Sg" \
    $FORCE_FLAG

# ── Stage 6c: dinucleotide proportions (2Sg) ─────────────────────────────────
log_info "Stage 6c: dinucleotide proportions — 2Sg"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SEL" \
    --samples "$ROOT_DIR/config/samples_2Sg.tsv" \
    --outdir  "$ROOT_DIR/results/dinuc_2Sg" \
    $FORCE_FLAG

# ── Stage 7: YR/YC enrichment ────────────────────────────────────────────────
log_info "Stage 7: YR/YC enrichment"
$RSCRIPT "$ROOT_DIR/scripts/04_dinuc/02_yr_enrichment.R" \
    --samples "$SAMPLES" \
    --outdir  "$ROOT_DIR/results/figures" \
    --params  "$PARAMS"

# ── Stages 8+9: TSS end-type matrix, clustering, annotation ──────────────────
log_info "Stages 8+9: TSS end-type matrix and clustering"
bash "$ROOT_DIR/scripts/05_tss_end_types/run.sh" $FORCE_FLAG

log_info "All stages complete."
