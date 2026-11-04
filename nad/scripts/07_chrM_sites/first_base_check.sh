#!/usr/bin/env bash
# NAD branch step 7, diagnostic: what starts at chrM:58010 (+)?
#
# Most NAD captureSeq S reads at the 21S rRNA promoter start 2 nt downstream of
# the final A of ATATAAGTA (58008), on a genomic T at 58010. Step 2 strips the
# whole leading G run, and a genomic G at +1 cannot be told apart from the TdT
# tail, so an RNA starting on the G at 58009 (GTAAAA...) maps to 58010. This
# script collects the evidence for each + strand 5' end in a window around
# 58008 and, as a baseline, around 73704 (an A start with no G after it):
#   - reads, reads with a mismatched first base (mm5), the first aligned base;
#   - the length of the stripped G run (from the read name <id>_g<N>_<UMI>):
#     if 58010 reads come from 58009, their G run is on average one longer
#     than at A starts, and reads with g0 (no G stripped) should be absent;
#   - the most frequent first 12 bases of the trimmed read;
#   - for one WT S run, the raw (untrimmed) reads behind a sample of read IDs;
#   - the 21S rRNA (Q0158) gene start in the Ensembl GTF.
#
# Run on vavilov through sbatch:
#   sbatch --parsable -p main -J nad-fb58010 -c 2 --mem=4G -t 1:00:00 \
#          -o ~/logs/%x-%j.out --export=ALL \
#          --wrap "bash nad/scripts/07_chrM_sites/first_base_check.sh"
#
# Outputs (nad/results/07_chrM_sites/first_base/):
#   per_position.tsv   run strain library locus pos reads mm5 g0 mean_glen first_base_counts
#   glen_hist.tsv      run locus pos glen reads
#   top_prefix.tsv     run locus pos prefix12 reads
#   raw_reads.tsv      run pos read_id raw_seq trimmed_seq   (WT S rep1 only)
#   q0158_gtf.txt      GTF lines of the 21S rRNA gene
#
# Usage: bash nad/scripts/07_chrM_sites/first_base_check.sh [--params config/params.yaml]
#        [--raw-run SRR11236402] [--n-raw 15]

set -euo pipefail
export LC_ALL=C   # decimal point in awk printf, plain byte sort
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
RAW_RUN="SRR11236402"
N_RAW=15
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)  PARAMS="$2"; shift 2 ;;
        --raw-run) RAW_RUN="$2"; shift 2 ;;
        --n-raw)   N_RAW="$2"; shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

BAM_DIR="$NAD_DIR/data/bam"
RAW_DIR="$NAD_DIR/data/nad_captureseq"
GTF="$NAD_DIR/data/annotation/Saccharomyces_cerevisiae.R64-1-1.113.gtf.gz"
OUT="$NAD_DIR/results/07_chrM_sites/first_base"
mkdir -p "$OUT" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/07_first_base_$(date +%Y%m%d_%H%M%S).log") 2>&1

SAMTOOLS="$(conda run -n "$(yaml_get conda_env_tools "$PARAMS")" which samtools)"
[[ -x "$SAMTOOLS" ]] || die "samtools not found in conda env $(yaml_get conda_env_tools "$PARAMS")"

# Loci: name, window start, window end (1-based, + strand 5' ends).
LOCI="21S 58004 58014
COX2 73700 73710"

# ── one pass per run over the two windows ────────────────────────────────────
printf "run\tstrain\tlibrary\tlocus\tpos\treads\tmm5\tg0\tmean_glen\tfirst_base_counts\n" > "$OUT/per_position.tsv"
printf "run\tlocus\tpos\tglen\treads\n"   > "$OUT/glen_hist.tsv"
printf "run\tlocus\tpos\tprefix12\treads\n" > "$OUT/top_prefix.tsv"
: > "$OUT/.ids_tmp"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run strain library rep _rest; do
    bam="$BAM_DIR/$run.bam"
    [[ -s "$bam" ]] || { log_warn "missing $bam"; continue; }
    while read -r locus from to; do
        # Region padded so that reads starting in the window are all fetched;
        # only + strand reads (flag 16 unset) whose POS lies in the window count.
        "$SAMTOOLS" view -F 16 -q 255 "$bam" "chrM:${from}-${to}" \
        | awk -F'\t' -v run="$run" -v strain="$strain" -v lib="$library" -v locus="$locus" \
              -v from="$from" -v to="$to" -v out="$OUT" -v raw_run="$RAW_RUN" '
            $4 < from || $4 > to { next }
            {
                pos = $4; seq = $10
                n = split($1, a, "_"); g = a[n-1]; sub(/^g/, "", g); g += 0
                reads[pos]++; glen_sum[pos] += g
                if (g == 0) g0[pos]++
                hist[pos "\t" g]++
                fb[pos "\t" substr(seq, 1, 1)]++
                pre[pos "\t" substr(seq, 1, 12)]++
                # MD tag: a 5'"'"' mismatch shows as MD:Z:0<base>...
                for (i = 12; i <= NF; i++) if ($i ~ /^MD:Z:0[ACGTN]/) mm[pos]++
                if (run == raw_run) print run "\t" pos "\t" a[1] "\t" seq >> (out "/.ids_tmp")
            }
            END {
                for (p = from; p <= to; p++) {
                    if (!(p in reads)) continue
                    fbs = ""
                    for (k in fb) { split(k, b, "\t"); if (b[1] == p) fbs = fbs b[2] ":" fb[k] " " }
                    printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%s\n", run, strain, lib, locus, p,
                           reads[p], mm[p] + 0, g0[p] + 0, glen_sum[p] / reads[p], fbs >> (out "/per_position.tsv")
                }
                for (k in hist) print run "\t" locus "\t" k "\t" hist[k] >> (out "/glen_hist.tsv")
                for (k in pre) if (pre[k] >= 3) print run "\t" locus "\t" k "\t" pre[k] >> (out "/top_prefix.tsv")
            }'
    done <<< "$LOCI"
    log_info "$run done"
done

# ── raw reads behind a sample of IDs (one WT S run) ──────────────────────────
# The read name in the BAM carries the original SRA id before the first "_"
# that step 2 added; SRA ids themselves contain a dot (SRRxxx.N), not "_".
printf "run\tpos\tread_id\traw_seq\ttrimmed_seq\n" > "$OUT/raw_reads.tsv"
raw="$RAW_DIR/$RAW_RUN.fastq.gz"
if [[ -s "$raw" ]]; then
    awk -F'\t' -v n="$N_RAW" '($2 == 58008 || $2 == 58010 || $2 == 58011) && c[$2]++ < n' \
        "$OUT/.ids_tmp" > "$OUT/.ids_pick"
    log_info "looking up $(wc -l < "$OUT/.ids_pick") reads in $raw"
    gzip -dc "$raw" | awk -F'\t' -v OFS='\t' '
        NR == FNR { want[$3] = $1 "\t" $2; trimmed[$3] = $4; next }
        FNR % 4 == 1 { split($0, h, " "); id = substr(h[1], 2); keep = (id in want); next }
        FNR % 4 == 2 && keep { print want[id], id, $0, trimmed[id] }
    ' "$OUT/.ids_pick" - | sort -t$'\t' -k2,2n >> "$OUT/raw_reads.tsv"
else
    log_warn "raw FASTQ not found: $raw (raw_reads.tsv stays empty)"
fi
rm -f "$OUT/.ids_tmp" "$OUT/.ids_pick"

# ── 21S rRNA annotation ──────────────────────────────────────────────────────
if [[ -s "$GTF" ]]; then
    gzip -dc "$GTF" | awk -F'\t' '$1 == "Mito" && $3 ~ /gene|exon|transcript/ && $9 ~ /Q0158|21S/' \
        | cut -f1,3,4,5,7,9 > "$OUT/q0158_gtf.txt"
else
    log_warn "GTF not found: $GTF"
fi

# ── summary to the log ───────────────────────────────────────────────────────
log_info "21S rRNA in the GTF:"; cut -f1-5 "$OUT/q0158_gtf.txt" 2>/dev/null || true
log_info "per position, pooled over runs by strain and library (reads, mm5, g0, mean G run):"
awk -F'\t' 'NR > 1 { k = $2 "\t" $3 "\t" $4 "\t" $5; r[k] += $6; m[k] += $7; z[k] += $8; g[k] += $9 * $6 }
            END { for (k in r) printf "%s\t%d\t%d\t%d\t%.2f\n", k, r[k], m[k], z[k], g[k] / r[k] }' \
    "$OUT/per_position.tsv" | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4n \
    | (printf "strain\tlibrary\tlocus\tpos\treads\tmm5\tg0\tmean_glen\n"; cat) | column -t
log_info "raw reads (WT S, sample):"
column -t "$OUT/raw_reads.tsv" | head -n $((3 * N_RAW + 1))
