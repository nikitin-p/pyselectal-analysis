#!/usr/bin/env bash
# Build a UCSC track hub (sacCer3) in results/tracks/, grouped by strand.
#
# The hub has two composite tracks, F (forward, + strand) and R (reverse, -
# strand), so that each group is selected as a whole. Each holds, as views:
#   - the eight STAR Signal.Unique bigWigs for that strand (four conditions x
#     two replicates), blue for F and red for R;
#   - the distclu TSS clusters of that strand coloured by hclust group, and the
#     same clusters coloured by dominant 5' end type, as bigBed 9 with itemRgb;
#   - the enhancers, which are unstranded and so appear in both groups;
#   - a per-base letter track at the top: the genome sequence in F and its
#     complement in R (scripts/99_additional/make_sequence_bed.py).
#
# STAR's str1 is the forward strand and str2 the reverse one for this library
# protocol. That was checked, not assumed: at the deepest + strand CTSS of
# yeast_cage_ana_r1 (chrII:477649) str1 carries the signal and str2 is empty,
# and at the deepest - strand CTSS (chrXII:370141) the reverse holds.
#
# Every file the hub names is placed next to hub.txt under the same basename
# it has on the cluster, so the folder can be uploaded to Hub Space or copied
# to a web server as it is. The TSS-cluster BEDs come from
# scripts/05_tss_end_types/11_export_cluster_bed.R.
#
# Usage (from project root):
#   bash scripts/99_additional/make_hub.sh [--bigwig_dir results/bigwig] \
#       [--outdir results/tracks] [--fai data/sacCer3.fa.fai] [--bedtobigbed PATH]

set -euo pipefail
export LC_ALL=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"

BW_DIR="$ROOT_DIR/results/bigwig"
OUTDIR="$ROOT_DIR/results/tracks"
FAI="$ROOT_DIR/data/sacCer3.fa.fai"
BEDTOBIGBED="$(command -v bedToBigBed || true)"
EMAIL="nikipavel14@gmail.com"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bigwig_dir)  BW_DIR="$2";      shift 2 ;;
        --outdir)      OUTDIR="$2";      shift 2 ;;
        --fai)         FAI="$2";         shift 2 ;;
        --bedtobigbed) BEDTOBIGBED="$2"; shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done
[[ -n "$BEDTOBIGBED" && -x "$BEDTOBIGBED" ]] || die "bedToBigBed not found (pass --bedtobigbed PATH)"
[[ -f "$FAI" ]] || die "Genome index not found: $FAI"

# Lookups are functions rather than associative arrays, because macOS ships
# bash 3.2, which has none.
CONDITIONS=(Ana Eth Glu Nit)
cond_long() { case "$1" in Ana) echo Anaerobic ;; Eth) echo Ethanol ;;
                           Glu) echo Glucose ;;   Nit) echo Nitrogen ;; esac; }
# The signal is coloured by strand, blue forward and red reverse, so the two
# groups read apart at a glance when both are on.
strand_color() { case "$1" in F) echo 33,102,172 ;; R) echo 178,24,43 ;; esac; }
str_of()  { case "$1" in F) echo 1 ;; R) echo 2 ;; esac; }
sign_of() { case "$1" in F) echo + ;; R) echo - ;; esac; }
name_of() { case "$1" in F) echo forward ;; R) echo reverse ;; esac; }
CLUSTER_SETS=(hclust_group dominant_type)
# UCSC truncates a shortLabel beyond 17 characters in the track list.
set_label() { case "$1" in hclust_group) echo "TSS group" ;;
                           dominant_type) echo "TSS dominant type" ;; esac; }
set_long() { case "$1" in
    hclust_group)  echo "distclu TSS clusters coloured by end-type profile group (Stage 9)" ;;
    dominant_type) echo "distclu TSS clusters coloured by dominant 5' end type at 0.7 (Stage 17); grey = not dominated" ;;
    esac; }

mkdir -p "$OUTDIR"

# ── files ─────────────────────────────────────────────────────────────────────
for c in "${CONDITIONS[@]}"; do
    for r in 1 2; do
        for s in 1 2; do
            f="$BW_DIR/${c}_Rep${r}.Signal.Unique.str${s}.out.wig.bw"
            [[ -f "$f" ]] || die "bigWig not found: $f"
            cp "$f" "$OUTDIR/"
        done
    done
done
log_info "copied 16 STAR bigWigs into $OUTDIR"

sizes="$(mktemp)"
trap 'rm -f "$sizes"' EXIT
cut -f1,2 "$FAI" > "$sizes"

for set in "${CLUSTER_SETS[@]}"; do
    bed="$OUTDIR/tss_clusters_by_${set}.bed"
    [[ -f "$bed" ]] || die "TSS-cluster BED not found: $bed (run 11_export_cluster_bed.R)"
    for g in F R; do
        body="$(mktemp)"
        awk -F'\t' -v s="$(sign_of $g)" '!/^track/ && $6 == s' "$bed" \
            | sort -k1,1 -k2,2n > "$body"
        "$BEDTOBIGBED" -type=bed9 "$body" "$sizes" "$OUTDIR/tss_clusters_by_${set}_${g}.bb"
        log_info "wrote tss_clusters_by_${set}_${g}.bb ($(wc -l < "$body" | tr -d ' ') clusters)"
        rm -f "$body"
    done
done

# Enhancers carry no strand, so both groups draw the one file. It is the
# UCSC-clean version written by ucsc_tracks.sh (integer scores, sorted).
enh_src="$OUTDIR/ucsc/bigbed/enhancers.bb"
[[ -f "$enh_src" ]] || die "enhancers bigBed not found: $enh_src (run ucsc_tracks.sh first)"
cp "$enh_src" "$OUTDIR/enhancers.bb"
log_info "copied enhancers.bb"

# Per-base letter tracks: the forward sequence heads F and its complement heads
# R, so each group can be read against the genome at base level. They are the
# largest files in the hub (one item per base) and depend only on the genome,
# so an existing pair is reused.
for g in F R; do
    bb="$OUTDIR/sequence_${g}.bb"
    if [[ -f "$bb" ]]; then
        log_info "sequence_${g}.bb exists, reusing it"
        continue
    fi
    body="$(mktemp)"
    python3 "$SCRIPT_DIR/make_sequence_bed.py" "${FAI%.fai}" "$body" \
        $([[ $g == R ]] && echo --complement)
    "$BEDTOBIGBED" -type=bed4 "$body" "$sizes" "$bb"
    rm -f "$body"
    log_info "wrote sequence_${g}.bb"
done

# ── hub.txt ───────────────────────────────────────────────────────────────────
hub="$OUTDIR/hub.txt"
{
cat <<EOF
#####
# hub config
#####

hub yeastCageHub
shortLabel Yeast CAGE 5prime
longLabel Yeast CAGE 5-end signal and TSS clusters across conditions (Ana, Eth, Glu, Nit), 2 replicates each, grouped by strand
descriptionUrl hubDescription.html
useOneFile on
email $EMAIL

#####
# genome config
#####

genome sacCer3
EOF

prio=1
for g in F R; do
    s="$(str_of $g)"
    # One composite per strand, so the whole group is switched on or off from a
    # single settings page (allButtonPair adds the select-all / clear-all pair).
    # Its members differ in type, which a composite allows only through views:
    # each view carries the type and display settings its subtracks share.
cat <<EOF

#####
# $g: $(name_of $g) ($(sign_of $g)) strand
#####

track $g
compositeTrack on
shortLabel $g
longLabel $(name_of $g | awk '{print toupper(substr($0,1,1)) substr($0,2)}') ($(sign_of $g)) strand: CAGE signal, TSS clusters and enhancers
subGroup1 view Views sequence=Sequence signal=Signal clusters=TSS_clusters enhancers=Enhancers
allButtonPair on
dragAndDrop subTracks
type bed 3
visibility full
priority $prio
html hubDescription.html

    track ${g}_view_sequence
    parent $g
    view sequence
    shortLabel $g sequence
    type bigBed 4
    color 255,255,255
    labelOnFeature on
    visibility dense

        track ${g}_sequence
        parent ${g}_view_sequence on
        subGroups view=sequence
        bigDataUrl sequence_${g}.bb
        shortLabel $([[ $g == F ]] && echo "Sequence +" || echo "Complement -")
        longLabel $([[ $g == F ]] && echo "Genome sequence, forward strand (one letter per base)" || echo "Complement of the genome sequence, in genome order; read right to left for the reverse strand 5' to 3'")
        type bigBed 4
        color 255,255,255
        labelOnFeature on
        priority 0

    track ${g}_view_signal
    parent $g
    view signal
    shortLabel $g signal
    type bigWig
    autoScale on
    visibility full

EOF
    p=1
    for c in "${CONDITIONS[@]}"; do
        for r in 1 2; do
cat <<EOF
        track ${g}_${c}_Rep${r}
        parent ${g}_view_signal on
        subGroups view=signal
        bigDataUrl ${c}_Rep${r}.Signal.Unique.str${s}.out.wig.bw
        shortLabel $c Rep$r $g
        longLabel $(cond_long $c) replicate $r, $(name_of $g) strand (STAR Signal.Unique.str$s)
        type bigWig
        color $(strand_color $g)
        priority $p

EOF
            p=$((p + 1))
        done
    done
cat <<EOF
    track ${g}_view_clusters
    parent $g
    view clusters
    shortLabel $g TSS clusters
    type bigBed 9
    itemRgb on
    visibility pack

EOF
    for set in "${CLUSTER_SETS[@]}"; do
cat <<EOF
        track ${g}_tss_${set}
        parent ${g}_view_clusters on
        subGroups view=clusters
        bigDataUrl tss_clusters_by_${set}_${g}.bb
        shortLabel $(set_label $set)
        longLabel $(set_long $set), $(name_of $g) strand
        type bigBed 9
        itemRgb on
        priority $p

EOF
        p=$((p + 1))
    done
cat <<EOF
    track ${g}_view_enhancers
    parent $g
    view enhancers
    shortLabel $g enhancers
    type bigBed 8
    visibility dense

        track ${g}_enhancers
        parent ${g}_view_enhancers on
        subGroups view=enhancers
        bigDataUrl enhancers.bb
        shortLabel Enhancers
        longLabel Enhancers (unstranded; the same regions are shown in both F and R)
        type bigBed 8
        priority $p
EOF
    prio=$((prio + 1))
done
} > "$hub"
log_info "wrote $hub"

# ── description page ──────────────────────────────────────────────────────────
# One page serves the hub and every track, since the tracks differ only in
# strand and condition; it carries the two colour legends, which a BED cannot.
legend_rows() {
    awk -F'\t' -v key="$2" -v lab="$3" 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
        { printf "<tr><td style=\"background:%s;width:2em\"></td><td>%s</td><td>%s</td><td align=\"right\">%s</td></tr>\n",
                 $h["colour"], $h[key], (lab == "" ? "" : $h[lab]), $h["n_clusters"] }' "$1"
}
strand_rows() {
    for g in F R; do
        printf '<tr><td style="background:rgb(%s);width:2em"></td><td>%s</td></tr>\n' "$(strand_color $g)" "$g, $(name_of $g) strand"
    done
}
page="$OUTDIR/hubDescription.html"
cat > "$page" <<HTML
<h2>Yeast CAGE 5′-end signal and TSS clusters</h2>
<p>This hub shows CAGE data from four growth conditions (Börlin et al. 2019), with two replicates each, on sacCer3. The tracks form two groups by strand: <b>F</b> holds the forward (+) strand and <b>R</b> holds the reverse (−) strand. Each group is one composite track, so its settings page switches all of its tracks on or off at once.</p>
<h3>Sequence</h3>
<p>The first track of each group spells out the genome one letter per base, in the style of the Base Position track. In F it is the forward sequence. In R it is the complement written in genome order, so reading it right to left gives the reverse strand 5′ to 3′, which is the strand the R signal and clusters come from. The letters appear once the view is zoomed to about base level; at wider views the track shows only its item density.</p>
<h3>Signal</h3>
<p>Each signal track is the STAR <code>Signal.Unique</code> coverage of uniquely mapped reads for one library and strand. STAR's <code>str1</code> is the forward strand and <code>str2</code> is the reverse strand for this protocol. The colour gives the strand, and the track name gives the condition and replicate.</p>
<table>
$(strand_rows)
</table>
<h3>TSS clusters</h3>
<p>The clusters are distclu tag clusters (maximum gap 20 bp, at least 5 pooled reads per position). Each box spans the cluster, the thick part marks its dominant TSS, and the item name ends with the initiator dinucleotide.</p>
<p><b>TSS group</b> colours each cluster by its end-type profile group from hierarchical clustering. Clusters present in fewer than two samples were not grouped and are absent from this track.</p>
<table>
<tr><th></th><th>Group</th><th>Mean profile</th><th>Clusters</th></tr>
$(legend_rows "$OUTDIR/tss_clusters_by_hclust_group_legend.tsv" hclust_group label)
</table>
<p><b>TSS dominant type</b> colours each cluster by the 5′ end type that dominates it at a threshold of 0.7, and clusters that no type dominates are grey. The item name gives the type and its share of the cluster's reads.</p>
<table>
<tr><th></th><th>Dominant type</th><th></th><th>Clusters</th></tr>
$(legend_rows "$OUTDIR/tss_clusters_by_dominant_type_legend.tsv" dom_label "")
</table>
<h3>Enhancers</h3>
<p>The enhancer regions carry no strand, so the same track appears in both F and R. Each item name gives the region followed by its original score.</p>
HTML
log_info "wrote $page"
