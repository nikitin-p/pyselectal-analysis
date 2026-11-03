#!/usr/bin/env bash
# NAD branch step 3: map trimmed NAD captureSeq reads to sacCer3 with STAR.
#
# The genome FASTA and GTF are the ones the CAGE pipeline used
# (genome_fasta / annotation_gtf in the params file), so that NAD 5' ends and
# CAGE CTSSs share coordinates exactly. The index is built once into
# nad/data/star_index_sacCer3 (yeast: --genomeSAindexNbases 10).
#
# Mapping differs from Stage 1 in one deliberate way: --alignEndsType
# Extend5pOfRead1 forbids a 5' soft clip, so the first aligned base is the
# RNA 5' end (the TdT G tail was already stripped in step 2). A non-templated
# 5' base therefore shows up as a mismatch, which step 4 records, instead of
# being hidden in a clip. Only unique alignments are kept (MAPQ 255).
#
# Outputs: nad/data/bam/<run>.bam (+ .bai), nad/results/03_map/<run>.Log.final.out,
#          nad/results/03_map/summary.tsv
# Idempotent: skips runs whose BAM exists unless --force.
#
# Usage: bash nad/scripts/03_map/run.sh [--params config/params.yaml]
#        [--manifest nad/scripts/00_download/priority1_runs.tsv] [--threads 10] [--force]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-10}"
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)   PARAMS="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

abs() { [[ "$1" == /* ]] && printf '%s' "$1" || printf '%s' "$ROOT_DIR/$1"; }
FASTA="$(abs "$(yaml_get genome_fasta "$PARAMS")")"
GTF="$(abs "$(yaml_get annotation_gtf "$PARAMS")")"
[[ -s "$FASTA" ]] || die "genome FASTA not found: $FASTA (copy data/sacCer3.fa from the Mac)"
[[ -s "$GTF"   ]] || die "annotation GTF not found: $GTF (copy data/sacCer3.ensGene.gtf from the Mac)"

TRIM_DIR="$NAD_DIR/data/trimmed"
IDX_DIR="$NAD_DIR/data/star_index_sacCer3"
BAM_DIR="$NAD_DIR/data/bam"
RES_DIR="$NAD_DIR/results/03_map"
mkdir -p "$BAM_DIR" "$RES_DIR" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/03_map_$(date +%Y%m%d_%H%M%S).log") 2>&1

STAR="$(conda run -n "$(yaml_get conda_env_star "$PARAMS")" which STAR)"
SAMTOOLS="$(conda run -n "$(yaml_get conda_env_tools "$PARAMS")" which samtools)"
log_info "STAR $("$STAR" --version), $("$SAMTOOLS" --version | head -1), threads $THREADS"

# ── index ─────────────────────────────────────────────────────────────────────
if [[ ! -s "$IDX_DIR/SA" ]]; then
    log_info "building STAR index → $IDX_DIR"
    mkdir -p "$IDX_DIR"
    "$STAR" --runMode genomeGenerate --runThreadN "$THREADS" \
        --genomeDir "$IDX_DIR" --genomeFastaFiles "$FASTA" \
        --sjdbGTFfile "$GTF" --sjdbOverhang 85 \
        --genomeSAindexNbases 10 \
        --outFileNamePrefix "$IDX_DIR/"
else
    log_info "index exists: $IDX_DIR"
fi

# ── map ───────────────────────────────────────────────────────────────────────
tail -n +2 "$MANIFEST" | cut -f1 | while read -r run; do
    fq="$TRIM_DIR/$run.fastq.gz"
    bam="$BAM_DIR/$run.bam"
    [[ -s "$fq" ]] || die "missing trimmed FASTQ: $fq (run step 02_trim first)"
    if [[ -s "$bam" && $FORCE -eq 0 ]]; then log_info "[skip] $run"; continue; fi

    tmp="$BAM_DIR/${run}_star_tmp"
    rm -rf "$tmp"; mkdir -p "$tmp"
    log_info "[map]  $run"
    # Short inserts (18 nt and up): cap mismatches at 10 % of the mapped length
    # instead of STAR's absolute 10; yeast introns are < 1.1 kb.
    "$STAR" --runThreadN "$THREADS" --genomeDir "$IDX_DIR" \
        --readFilesIn "$fq" --readFilesCommand zcat \
        --alignEndsType Extend5pOfRead1 \
        --outFilterMultimapNmax 1 --outSAMmultNmax 1 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMax 1100 \
        --outSAMattributes NH HI AS NM MD \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 3000000000 \
        --outFileNamePrefix "$tmp/" > /dev/null

    "$SAMTOOLS" view -b -q 255 -@ "$THREADS" "$tmp/Aligned.sortedByCoord.out.bam" -o "$bam.tmp"
    mv "$bam.tmp" "$bam"
    "$SAMTOOLS" index "$bam"
    cp "$tmp/Log.final.out" "$RES_DIR/$run.Log.final.out"
    rm -rf "$tmp"
    log_info "[ok]   $run: $(grep 'Uniquely mapped reads %' "$RES_DIR/$run.Log.final.out" | tr -s ' \t' ' ')"
done

# ── summary ───────────────────────────────────────────────────────────────────
{
    printf "run\tinput_reads\tunique\tunique_pct\tmulti_pct\ttoo_many_loci_pct\tmismatch_pct\ttoo_short_pct\tother_pct\n"
    tail -n +2 "$MANIFEST" | cut -f1 | while read -r run; do
        f="$RES_DIR/$run.Log.final.out"; [[ -s "$f" ]] || continue
        g() { grep "$1" "$f" | head -1 | awk -F'|' '{gsub(/[ \t]/,"",$2); print $2}'; }
        # With --outFilterMultimapNmax 1 every multimapper lands in "too many
        # loci", so multi_pct is always 0 and the loss shows up in that column.
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run" \
            "$(g 'Number of input reads')" "$(g 'Uniquely mapped reads number')" \
            "$(g 'Uniquely mapped reads %')" "$(g '% of reads mapped to multiple loci')" \
            "$(g '% of reads mapped to too many loci')" \
            "$(g '% of reads unmapped: too many mismatches')" \
            "$(g '% of reads unmapped: too short')" "$(g '% of reads unmapped: other')"
    done
} > "$RES_DIR/summary.tsv"
column -t "$RES_DIR/summary.tsv"
