#!/usr/bin/env bash
# Stage 22: where the long 5' soft clips come from.
#
# A long clip is one longer than dom_max_clip (the "longclip" fine type). The
# stage isolates the clipped segment of every such read in the subsampled BAMs,
# places it in the genome (the end of exon 1 for a spliced read, or an exact
# match upstream on the read strand), runs the reversed clip through the same
# steps as a chance-match control, and annotates each read's inferred 5' end
# with its TSS cluster class and the 1Sg reads beside it.
#
# Steps:
#   1. Classify each clip's origin and annotate the inferred 5' end
#   2. Tables for each question and the figure 26_longclip_origins.pdf
#   3. Browser tracks: a BED12 (bigBed with --bedtobigbed) drawing each clip as a
#      thick block at its origin, and bigWigs per origin class
#   4. Match/mismatch BAMs: every long-clip read realigned with its clip at its
#      origin and an =/X CIGAR, a case-encoded read in the ms:Z tag (UPPER match,
#      lower mismatch), and a segment-split BAM coloured by YC (pysam)
#   5. M counterparts: unclipped reads at each long-clip 5' end, and clipped vs
#      split reads by exon-1 piece length (pysam)
#   6. Figure 27_longclip_m_counterparts.pdf
#   7. Location in gene structure, clip length per location, and 5' structure
#      (including the untemplated junk) against location: tables
#      lc_location_*.tsv and figure 28_longclip_location.pdf
#
# Usage:
#   bash scripts/14_longclip/run.sh [--params ...] [--samples ...] [--force] \
#       [--bedtobigbed PATH] [--conda-env NAME]
#
# --conda-env names the env that has pysam for step 4; the default is
# conda_env_pyselectal from the params file.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
SAMPLES="$ROOT_DIR/config/samples.tsv"
FORCE_FLAG=""
PY_ENV=""
BEDTOBIGBED="$(command -v bedToBigBed || true)"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)      PARAMS="$2";      shift 2 ;;
        --samples)     SAMPLES="$2";     shift 2 ;;
        --force)       FORCE_FLAG="--force"; shift ;;
        --bedtobigbed) BEDTOBIGBED="$2"; shift 2 ;;
        --conda-env)   PY_ENV="$2";      shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done
RESULTS_ROOT="$(resolve_results_root "$PARAMS" "$ROOT_DIR")"

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
ANNOTATION_GTF="$(resolve_path "$(yaml_get annotation_gtf "$PARAMS")")"
DOM_DIR="$ROOT_DIR/$(yaml_get results_dominance "$PARAMS")"
EXO_DIR="$ROOT_DIR/$(yaml_get results_exotic_origins "$PARAMS")"
LC_DIR="$ROOT_DIR/$(yaml_get results_longclip "$PARAMS")"
FIG_DIR="$RESULTS_ROOT/figures"

LOG="$ROOT_DIR/logs/22_longclip_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$(dirname "$LOG")"

log_info "Stage 22: origins of the long 5' soft clips"

log_info "  step 1: classify clip origins"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/01_classify.R" \
    --samples   "$SAMPLES" \
    --fasta     "$GENOME_FASTA" \
    --gtf       "$ANNOTATION_GTF" \
    --typed     "$RESULTS_ROOT/tss/typed_ctss.tsv.gz" \
    --matrix    "$RESULTS_ROOT/tss/tss_matrix.tsv" \
    --promoters "$RESULTS_ROOT/tss/promoter_clusters.tsv" \
    --fine      "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --dominance "$EXO_DIR/position_dominance.tsv.gz" \
    --outdir    "$LC_DIR" \
    --params    "$PARAMS" $FORCE_FLAG

log_info "  step 2: tables and figure"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/02_analysis.R" \
    --reads     "$LC_DIR/longclip_reads.tsv.gz" \
    --fine      "$DOM_DIR/fine_typed_ctss.tsv.gz" \
    --promoters "$RESULTS_ROOT/tss/promoter_clusters.tsv" \
    --fasta     "$GENOME_FASTA" \
    --gtf       "$ANNOTATION_GTF" \
    --outdir    "$LC_DIR" \
    --figdir    "$FIG_DIR" \
    --params    "$PARAMS" $FORCE_FLAG

log_info "  step 3: browser tracks"
[[ -z "$BEDTOBIGBED" ]] && log_warn "bedToBigBed not found; only the BED12 is written (pass --bedtobigbed PATH)"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/03_tracks.R" \
    --reads  "$LC_DIR/longclip_reads.tsv.gz" \
    --fai    "$GENOME_FASTA.fai" \
    --outdir "$LC_DIR/tracks" \
    ${BEDTOBIGBED:+--bedtobigbed "$BEDTOBIGBED"} $FORCE_FLAG

log_info "  step 4: match/mismatch BAMs"
[[ -z "$PY_ENV" ]] && PY_ENV="$(yaml_get conda_env_pyselectal "$PARAMS")"
log_to_file "$LOG" conda run -n "$PY_ENV" python "$SCRIPT_DIR/04_case_bams.py" \
    --samples  "$SAMPLES" \
    --reads    "$LC_DIR/longclip_reads.tsv.gz" \
    --fasta    "$GENOME_FASTA" \
    --helpers  "$SCRIPT_DIR/lc_helpers.R" \
    --max-clip "$(yaml_get dom_max_clip "$PARAMS" | tr -d '[:space:]')" \
    --outdir   "$LC_DIR/bam" $FORCE_FLAG

log_info "  step 5: M counterparts"
log_to_file "$LOG" conda run -n "$PY_ENV" python "$SCRIPT_DIR/05_m_counterparts.py" \
    --samples   "$SAMPLES" \
    --reads     "$LC_DIR/longclip_reads.tsv.gz" \
    --realigned "$LC_DIR/bam/longclip_realigned.bam" \
    --max-clip  "$(yaml_get dom_max_clip "$PARAMS" | tr -d '[:space:]')" \
    --outdir    "$LC_DIR" $FORCE_FLAG

log_info "  step 6: M-counterpart figure"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/06_m_counterparts_plot.R" \
    --indir  "$LC_DIR" \
    --genes  "$ROOT_DIR/$(yaml_get annotation_bed_yeast "$PARAMS")" \
    --figdir "$FIG_DIR" $FORCE_FLAG

log_info "  step 7: location summary"
log_to_file "$LOG" "$RSCRIPT" "$SCRIPT_DIR/07_location.R" \
    --reads     "$LC_DIR/longclip_reads.tsv.gz" \
    --gtf       "$ANNOTATION_GTF" \
    --promoters "$RESULTS_ROOT/tss/promoter_clusters.tsv" \
    --outdir    "$LC_DIR" \
    --figdir    "$FIG_DIR" \
    --params    "$PARAMS" $FORCE_FLAG

log_info "Stage 22 complete."
log_info "  tables → $LC_DIR"
log_info "  tracks → $LC_DIR/tracks"
log_info "  BAMs   → $LC_DIR/bam"
log_info "  figure → $FIG_DIR/26_longclip_origins.pdf, 27_longclip_m_counterparts.pdf, 28_longclip_location.pdf"
