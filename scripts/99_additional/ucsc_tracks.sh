#!/usr/bin/env bash
# Make the BED tracks in results/tracks/ loadable as UCSC custom tracks (sacCer3).
#
# The CAGEr BEDs extracted from results/yeast_tracks.tar.gz break three of
# UCSC's BED rules, which IGV tolerates:
#   - the score column holds a TPM, which is fractional and often above 1000,
#     whereas UCSC requires an integer in 0..1000;
#   - consensusClusters.bed carries a 7th column that is always 0, which UCSC
#     reads as a thickStart lying outside the feature, and its track line was
#     written with tabs instead of spaces;
#   - enhancers.bed has the same fractional score.
# Each file is rewritten into <outdir>: the score is rounded and capped at
# 1000, and the exact value moves into the name column (tpm=… for the CAGEr
# clusters, appended after "|" for the enhancers) so that nothing is lost. The
# two TSS-cluster tracks written by
# scripts/05_tss_end_types/11_export_cluster_bed.R already conform and are
# copied unchanged. All files are also concatenated into one upload file, since
# UCSC accepts several track lines in a single custom-track upload.
#
# Those .bed files load through the Custom Tracks page (hgCustom) only. Hub
# Space (My Data -> Track Hubs -> Hub Upload) rejects plain BED with "File type
# not supported" and accepts indexed binaries such as bigBed and bigWig, so when
# UCSC's bedToBigBed is available each track is also converted into
# <outdir>/bigbed/<stem>.bb, with its track line stripped and its rows sorted as
# bedToBigBed requires.
#
# Usage (from project root):
#   bash scripts/99_additional/ucsc_tracks.sh [--indir results/tracks] \
#       [--outdir results/tracks/ucsc] [--fai data/sacCer3.fa.fai]    \
#       [--bedtobigbed PATH]
#
# --bedtobigbed defaults to bedToBigBed on PATH; without it the bigBed step is
# skipped with a warning. The binary is at
# https://hgdownload.soe.ucsc.edu/admin/exe/ (macOSX.arm64 on Apple silicon).

set -euo pipefail
# awk honours LC_NUMERIC: under a comma-decimal locale (ru_RU here) it reads
# "17.806" as 17 and prints "17,000", which would silently truncate every score.
export LC_ALL=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"

INDIR="$ROOT_DIR/results/tracks"
OUTDIR=""
FAI="$ROOT_DIR/data/sacCer3.fa.fai"
BEDTOBIGBED="$(command -v bedToBigBed || true)"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --indir)  INDIR="$2";  shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --fai)    FAI="$2";    shift 2 ;;
        --bedtobigbed) BEDTOBIGBED="$2"; shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done
[[ -z "$OUTDIR" ]] && OUTDIR="$INDIR/ucsc"
mkdir -p "$OUTDIR"

# Round a score into UCSC's 0..1000 integer range.
AWK_SCORE='function ucsc_score(x) { x = int(x + 0.5); return (x > 1000) ? 1000 : ((x < 0) ? 0 : x) }'

# CAGEr tag clusters (BED12) and consensus clusters (BED7 -> BED6).
for f in "$INDIR"/*_tagClusters.bed "$INDIR"/consensusClusters.bed; do
    [[ -f "$f" ]] || continue
    stem="$(basename "$f" .bed)"
    awk -F'\t' -v OFS='\t' -v stem="$stem" "$AWK_SCORE"'
        /^track/ { printf "track name=\"%s\" description=\"CAGEr %s (score = TPM, capped at 1000; exact TPM in name)\" visibility=dense\n", stem, stem; next }
        {
            $4 = sprintf("tpm=%.3f", $5)
            $5 = ucsc_score($5)
            if (NF == 7) NF = 6
            print
        }' "$f" > "$OUTDIR/$stem.bed"
    log_info "wrote $OUTDIR/$stem.bed"
done

# Enhancers (BED8, no track line).
if [[ -f "$INDIR/enhancers.bed" ]]; then
    awk -F'\t' -v OFS='\t' "$AWK_SCORE"'
        BEGIN { print "track name=\"enhancers\" description=\"Enhancers (score rounded, capped at 1000; exact score in name)\" visibility=dense" }
        !/^track/ { $4 = sprintf("%s|%.3f", $4, $5); $5 = ucsc_score($5); print }
        ' "$INDIR/enhancers.bed" > "$OUTDIR/enhancers.bed"
    log_info "wrote $OUTDIR/enhancers.bed"
fi

# The pipeline's own TSS-cluster tracks already conform.
for f in "$INDIR"/tss_clusters_by_*.bed; do
    [[ -f "$f" ]] || continue
    cp "$f" "$OUTDIR/"
    log_info "copied $OUTDIR/$(basename "$f")"
done

# One file for a single upload, the pipeline's tracks first.
combined="$OUTDIR/all_tracks_sacCer3.bed"
cat "$OUTDIR"/tss_clusters_by_*.bed "$OUTDIR"/consensusClusters.bed \
    "$OUTDIR"/*_tagClusters.bed "$OUTDIR"/enhancers.bed > "$combined"
log_info "wrote $combined"

# ── bigBed for Hub Space ─────────────────────────────────────────────────────
if [[ -z "$BEDTOBIGBED" || ! -x "$BEDTOBIGBED" ]]; then
    log_warn "bedToBigBed not found; skipping the bigBed conversion (pass --bedtobigbed PATH)"
    exit 0
fi
[[ -f "$FAI" ]] || die "Genome index not found: $FAI"

BBDIR="$OUTDIR/bigbed"
mkdir -p "$BBDIR"
sizes="$BBDIR/sacCer3.chrom.sizes"
cut -f1,2 "$FAI" > "$sizes"

for f in "$OUTDIR"/*.bed; do
    stem="$(basename "$f" .bed)"
    [[ "$stem" == all_tracks_sacCer3 ]] && continue
    body="$BBDIR/$stem.sorted.bed"
    grep -v '^track' "$f" | sort -k1,1 -k2,2n > "$body"
    ncol=$(awk -F'\t' 'NR == 1 { print NF; exit }' "$body")
    "$BEDTOBIGBED" -type="bed$ncol" "$body" "$sizes" "$BBDIR/$stem.bb"
    rm -f "$body"
    log_info "wrote $BBDIR/$stem.bb (bed$ncol)"
done
