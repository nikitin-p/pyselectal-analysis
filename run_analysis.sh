#!/usr/bin/env bash
# Top-level runner: stages 1–9 (mapping → dinucleotides → TSS matrix).
#
# Usage:
#   bash run_analysis.sh [--force] [--skip-mapping] [--params FILE] [--samples FILE]
#
# Options:
#   --force         Re-run stages even if outputs exist
#   --skip-mapping  Skip stage 1 (STAR mapping); useful when BAMs already exist
#   --params        Parameter file           (default config/params.yaml)
#   --samples       Sample sheet             (default config/samples.tsv)
#
# Every output path is derived from the params file's results_root, so a second
# dataset is run by passing its own params and samples files rather than by
# editing this script. For the human LUHMES data:
#
#   bash run_analysis.sh --params config/params_human.yaml \
#                        --samples config/samples_human.tsv
#
# Prerequisites:
#   - FASTQ paths filled in the sample sheet
#   - STAR index built and named in the params file

set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
SKIP_MAPPING=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)        FORCE_FLAG="--force"; shift ;;
        --skip-mapping) SKIP_MAPPING=true;    shift ;;
        --params)       PARAMS="$2";  shift 2 ;;
        --samples)      SAMPLES="$2"; shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# Resolve every output location from the params file, so that two datasets never
# write into the same tree. The per-type sheets that stage 6a generates are
# stemmed on the sample sheet's own basename for the same reason.
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"
SHEET_PREFIX="$(basename "$SAMPLES" .tsv)"

CONDA_R="$(yaml_get conda_env_r "$PARAMS")"
RSCRIPT=$(conda run -n "$CONDA_R" which Rscript)

BAM_SUB="$RESULTS_ROOT/bam_subsampled"
BAM_SEL="$RESULTS_ROOT/bam_selected"

# ── Stage 1: STAR mapping ─────────────────────────────────────────────────────
if [[ "$SKIP_MAPPING" == true ]]; then
    log_info "Stage 1: STAR mapping — SKIPPED (--skip-mapping)"
else
    log_info "Stage 1: STAR mapping"
    bash "$ROOT_DIR/scripts/01_map/01_map.sh" \
        --samples "$SAMPLES" \
        --params  "$PARAMS" \
        $FORCE_FLAG
fi

# ── Stage 2: subsampling ──────────────────────────────────────────────────────
log_info "Stage 2: subsampling to equal depth"
bash "$ROOT_DIR/scripts/02_subsample/01_subsample.sh" --params "$PARAMS" $FORCE_FLAG

# ── Stage 3: pyselectal --count ───────────────────────────────────────────────
log_info "Stage 3: 5' end type counts"
bash "$ROOT_DIR/scripts/03_five_prime/01_count.sh" \
    --samples "$SAMPLES" \
    --params  "$PARAMS" \
    --outdir  "$RESULTS_ROOT/five_prime" \
    $FORCE_FLAG

# ── Stage 4: per-sample histograms + heatmap + chi-square ────────────────────
log_info "Stage 4: replicate heatmap and chi-square"
$RSCRIPT "$ROOT_DIR/scripts/03_five_prime/02_plot.R" \
    --counts_dir "$RESULTS_ROOT/five_prime" \
    --samples    "$SAMPLES" \
    --outdir     "$RESULTS_ROOT/figures" \
    --params     "$PARAMS"

# ── Stage 5: soft-clip nucleotide composition ─────────────────────────────────
log_info "Stage 5: soft-clip nucleotide composition"
$RSCRIPT "$ROOT_DIR/scripts/03_five_prime/03_softclip_composition.R" \
    --counts_dir "$RESULTS_ROOT/five_prime" \
    --samples    "$SAMPLES" \
    --outdir     "$RESULTS_ROOT/figures" \
    --params     "$PARAMS"

# ── Stage 6a: select per-type BAMs ───────────────────────────────────────────
log_info "Stage 6a: select 1Sg, 2Sgg, and 3Sggg BAMs"
bash "$ROOT_DIR/scripts/04_dinuc/00_select_bam.sh" \
    --bam_dir "$BAM_SUB" \
    --samples "$SAMPLES" \
    --params  "$PARAMS" \
    --sheet_prefix "$SHEET_PREFIX" \
    $FORCE_FLAG

# ── Stage 6b: dinucleotide proportions (1Sg) ─────────────────────────────────
log_info "Stage 6b: dinucleotide proportions — 1Sg"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SEL" \
    --samples "$ROOT_DIR/config/${SHEET_PREFIX}_1Sg.tsv" \
    --params  "$PARAMS" \
    --outdir  "$RESULTS_ROOT/dinuc_1Sg" \
    $FORCE_FLAG

# ── Stage 6c: dinucleotide proportions (2Sgg) ────────────────────────────────
log_info "Stage 6c: dinucleotide proportions — 2Sgg"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SEL" \
    --samples "$ROOT_DIR/config/${SHEET_PREFIX}_2Sgg.tsv" \
    --params  "$PARAMS" \
    --outdir  "$RESULTS_ROOT/dinuc_2Sgg" \
    $FORCE_FLAG

# ── Stage 6d: dinucleotide proportions (3Sggg) ───────────────────────────────
log_info "Stage 6d: dinucleotide proportions — 3Sggg"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SEL" \
    --samples "$ROOT_DIR/config/${SHEET_PREFIX}_3Sggg.tsv" \
    --params  "$PARAMS" \
    --outdir  "$RESULTS_ROOT/dinuc_3Sggg" \
    $FORCE_FLAG

# ── Stage 6e: dinucleotide proportions (all reads — unfiltered subsampled BAMs)
log_info "Stage 6e: dinucleotide proportions — all reads"
bash "$ROOT_DIR/scripts/04_dinuc/run.sh" \
    --bam_dir "$BAM_SUB" \
    --samples "$SAMPLES" \
    --params  "$PARAMS" \
    --outdir  "$RESULTS_ROOT/dinuc_all" \
    $FORCE_FLAG

# ── Stage 7: YR/YC enrichment ────────────────────────────────────────────────
log_info "Stage 7: YR/YC enrichment"
$RSCRIPT "$ROOT_DIR/scripts/04_dinuc/02_yr_enrichment.R" \
    --dinuc_1sg   "$RESULTS_ROOT/dinuc_1Sg/all_dinuc_proportions.tsv" \
    --dinuc_2sgg  "$RESULTS_ROOT/dinuc_2Sgg/all_dinuc_proportions.tsv" \
    --dinuc_3sggg "$RESULTS_ROOT/dinuc_3Sggg/all_dinuc_proportions.tsv" \
    --dinuc_all   "$RESULTS_ROOT/dinuc_all/all_dinuc_proportions.tsv" \
    --samples "$SAMPLES" \
    --outdir  "$RESULTS_ROOT/figures" \
    --params  "$PARAMS"

# ── Stages 8+9: TSS end-type matrix, clustering, annotation ──────────────────
log_info "Stages 8+9: TSS end-type matrix and clustering"
bash "$ROOT_DIR/scripts/05_tss_end_types/run.sh" --params "$PARAMS" --samples "$SAMPLES" $FORCE_FLAG

# ── Stage 17: single-type-dominated TSSs ─────────────────────────────────────
log_info "Stage 17: single-type-dominated TSSs"
bash "$ROOT_DIR/scripts/10_dominant_endtype/run.sh" --params "$PARAMS" --samples "$SAMPLES" $FORCE_FLAG

# ── Stage 18: A-8 saturation of <5' end type, dinucleotide> pairs ────────────
log_info "Stage 18: A-8 saturation of <5' end type, dinucleotide> pairs"
bash "$ROOT_DIR/scripts/11_a8_saturation/run.sh" --params "$PARAMS" --samples "$SAMPLES" $FORCE_FLAG

# ── Stage 19: origins of the exotic 5' end types ─────────────────────────────
# Reads Stage 17's tables, so it must follow it. It takes no sample sheet,
# because every input is already pooled over the libraries.
log_info "Stage 19: origins of the exotic 5' end types"
bash "$ROOT_DIR/scripts/12_exotic_origins/run.sh" --params "$PARAMS" $FORCE_FLAG

# ── Stage 14: HTML analysis report ───────────────────────────────────────────
# The report's prose, its figures and its synthesis are written about the yeast
# run specifically, so it is rendered only for the default results tree. Running
# it against another dataset would produce yeast conclusions illustrated with
# that dataset's numbers, which is worse than producing no report at all.
if [[ "$RESULTS_ROOT" == "$ROOT_DIR/results" ]]; then
    log_info "Stage 14: HTML analysis report"
    bash "$ROOT_DIR/scripts/07_report/run.sh" $FORCE_FLAG
else
    log_info "Stage 14: HTML analysis report — SKIPPED (report is specific to the yeast run)"
fi

log_info "All stages complete."
