#!/usr/bin/env bash
# NAD branch step 7, diagnostic: are NAD captureSeq reads from chrM 6470 and
# 46079 lost between the raw FASTQ and the BAM?
#
# In CAGE, 6470 (15S promoter) and 46079 (ATP9) carry 30-37 % 1Sa, yet the
# step-3 BAMs hold at most three S reads per replicate there. A genome-only
# check found that the 18-mer starting at +1 of both sites occurs once in
# sacCer3 on either strand, so an error-free read maps uniquely. This script
# looks for the reads themselves at three points of the chain, per run and site:
#   bam      reads in the step-3 BAM (MAPQ 255) whose 5' end is +1 on the site strand;
#   trimmed  reads in the step-2 FASTQ whose bases 2..18 equal the genome +2..+18
#            (any first base, so a mismatched 5' base is caught), with their first
#            base and median length;
#   raw      reads in the downloaded FASTQ that contain the genomic 18-mer +1..+18
#            anywhere, with the prefix before it (expected: the TdT G run).
# If trimmed or raw is much larger than bam, reads are lost in mapping or
# trimming; if all three are near zero, the RNA is absent from the library.
# Controls: 58008, 73704, 13357, 32217 (detected, unique), 6548 (not detected)
# and 54840(-) (detected, but its 18-mer occurs 7 times in the genome, 82329
# among them, so its trimmed and raw counts pool all seven loci while bam counts
# the one position; 82329 is left out because it shares the 18-mer).
#
# Run on vavilov through sbatch:
#   sbatch --parsable -p main -J nad-rawchk -c 2 --mem=4G -t 4:00:00 \
#          -o ~/logs/%x-%j.out --export=ALL \
#          --wrap "bash nad/scripts/07_chrM_sites/raw_read_check.sh"
#
# Outputs (nad/results/07_chrM_sites/raw_read_check/):
#   sites.tsv        site strand seq18
#   counts.tsv       run strain library site bam trimmed trimmed_A trimmed_mm5 trimmed_median_len raw raw_allG_prefix
#   raw_prefix.tsv   run site prefix reads    (prefix before the 18-mer in raw reads, top entries)
#
# Usage: bash nad/scripts/07_chrM_sites/raw_read_check.sh [--params config/params.yaml]
#        [--manifest ...] [--runs SRR1,SRR2] [--no-bam]

set -euo pipefail
export LC_ALL=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAD_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
ROOT_DIR="$(cd "$NAD_DIR/.." && pwd)"
source "$ROOT_DIR/scripts/utils/log.sh"
source "$ROOT_DIR/scripts/utils/parse_yaml.sh"

PARAMS="$ROOT_DIR/config/params.yaml"
MANIFEST="$NAD_DIR/scripts/00_download/priority1_runs.tsv"
RUNS=""; NO_BAM=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)   PARAMS="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --runs)     RUNS="$2"; shift 2 ;;
        --no-bam)   NO_BAM=1; shift ;;
        *) die "Unknown argument: $1" ;;
    esac
done

abs() { [[ "$1" == /* ]] && printf '%s' "$1" || printf '%s' "$ROOT_DIR/$1"; }
FASTA="$(abs "$(yaml_get genome_fasta "$PARAMS")")"
[[ -s "$FASTA" ]] || die "genome FASTA not found: $FASTA"
BAM_DIR="$NAD_DIR/data/bam"
TRIM_DIR="$NAD_DIR/data/trimmed"
RAW_DIR="$NAD_DIR/data/nad_captureseq"
OUT="$NAD_DIR/results/07_chrM_sites/raw_read_check"
mkdir -p "$OUT" "$NAD_DIR/logs"
exec > >(tee -a "$NAD_DIR/logs/07_raw_read_check_$(date +%Y%m%d_%H%M%S).log") 2>&1

SAMTOOLS=""
if [[ $NO_BAM -eq 0 ]]; then
    SAMTOOLS="$(conda run -n "$(yaml_get conda_env_tools "$PARAMS")" which samtools)"
    [[ -x "$SAMTOOLS" ]] || die "samtools not found in conda env $(yaml_get conda_env_tools "$PARAMS")"
fi

# ── site 18-mers from the genome (strand-aware, FASTA is soft-masked) ───────────
SITES="6470 +
46079 +
6548 +
58008 +
73704 +
13357 +
32217 +
54840 -"
printf '%s\n' "$SITES" | awk -v fa="$FASTA" '
    BEGIN { while ((getline l < fa) > 0) {
                if (l ~ /^>/) { split(substr(l, 2), h, /[ \t]/); on = (h[1] == "chrM"); continue }
                if (on) m = m toupper(l) }
            if (m == "") { print "chrM not found in " fa > "/dev/stderr"; exit 1 }
            c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A"; c["N"] = "N" }
    { p = $1; s = $2
      if (s == "+") q = substr(m, p, 18)
      else { f = substr(m, p - 17, 18); q = ""; for (i = 18; i >= 1; i--) q = q c[substr(f, i, 1)] }
      print p "\t" s "\t" q }' > "$OUT/sites.tsv"
log_info "sites (18-mer from +1, RNA sense):"; column -t "$OUT/sites.tsv"

# ── per run ─────────────────────────────────────────────────────────────────────
printf "run\tstrain\tlibrary\tsite\tbam\ttrimmed\ttrimmed_A\ttrimmed_mm5\ttrimmed_median_len\traw\traw_allG_prefix\n" > "$OUT/counts.tsv"
printf "run\tsite\tprefix\treads\n" > "$OUT/raw_prefix.tsv"

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r run strain library _rest; do
    if [[ -n "$RUNS" && ",$RUNS," != *",$run,"* ]]; then continue; fi
    trim="$TRIM_DIR/$run.fastq.gz"; raw="$RAW_DIR/$run.fastq.gz"; bam="$BAM_DIR/$run.bam"
    [[ -s "$trim" && -s "$raw" ]] || { log_warn "missing FASTQ for $run, skipped"; continue; }
    log_info "[run] $run ($strain $library)"

    # BAM: 5' end of + reads is POS; of - reads, POS + reference span - 1.
    : > "$OUT/.bam_tmp"
    if [[ $NO_BAM -eq 0 ]]; then
        [[ -s "$bam" ]] || die "missing BAM: $bam"
        "$SAMTOOLS" view -q 255 "$bam" chrM | awk -F'\t' -v sf="$OUT/sites.tsv" '
                BEGIN { while ((getline l < sf) > 0) { split(l, a, "\t"); want[a[1] "\t" a[2]] = a[1] } }
                { minus = int($2 / 16) % 2
                  if (!minus) { e = $4; st = "+" }
                  else { span = 0; cg = $6
                         while (match(cg, /^[0-9]+[MIDNSHP=X]/)) {
                             len = substr(cg, 1, RLENGTH - 1) + 0; op = substr(cg, RLENGTH, 1)
                             if (op ~ /[MDN=X]/) span += len
                             cg = substr(cg, RLENGTH + 1) }
                         e = $4 + span - 1; st = "-" }
                  k = e "\t" st; if (k in want) n[want[k]]++ }
                END { for (k in want) print want[k] "\t" n[want[k]] + 0 }' > "$OUT/.bam_tmp"
    fi

    # Trimmed FASTQ: key = read bases 2..18 against genome +2..+18.
    trim_res="$(gzip -dc "$trim" | awk -v sf="$OUT/sites.tsv" '
        BEGIN { while ((getline l < sf) > 0) { split(l, a, "\t"); key[substr(a[3], 2, 17)] = a[1]; first[a[1]] = substr(a[3], 1, 1) } }
        NR % 4 == 2 && length($0) >= 18 {
            k = substr($0, 2, 17)
            if (k in key) { s = key[k]; n[s]++; b = substr($0, 1, 1)
                            if (b == "A") nA[s]++
                            if (b != first[s]) mm[s]++
                            lens[s] = lens[s] " " length($0) } }
        END { for (k in key) { s = key[k]; med = "NA"
                  if (n[s] > 0) { m = split(substr(lens[s], 2), L, " ")
                      for (i = 2; i <= m; i++) { v = L[i]; j = i - 1; while (j >= 1 && L[j] > v) { L[j + 1] = L[j]; j-- }; L[j + 1] = v }
                      med = L[int((m + 1) / 2)] }
                  print s "\t" n[s] + 0 "\t" nA[s] + 0 "\t" mm[s] + 0 "\t" med } }')"

    # Raw FASTQ: 18-mer anywhere; prefix before it.
    raw_res="$(gzip -dc "$raw" | awk -v sf="$OUT/sites.tsv" -v run="$run" -v pf="$OUT/raw_prefix.tsv" '
        BEGIN { while ((getline l < sf) > 0) { split(l, a, "\t"); q[a[1]] = a[3] } }
        NR % 4 == 2 { for (s in q) { i = index($0, q[s]); if (i > 0) {
                          n[s]++; pre = substr($0, 1, i - 1)
                          if (pre ~ /^G*$/) g[s]++
                          if (length(pre) > 12) pre = substr(pre, 1, 12) "+"
                          P[s "\t" pre]++ } } }
        END { for (s in q) print s "\t" n[s] + 0 "\t" g[s] + 0
              for (k in P) if (P[k] >= 2) print run "\t" k "\t" P[k] >> pf }')"

    while IFS=$'\t' read -r site _strand _seq; do
        t="$(awk -F'\t' -v s="$site" '$1 == s { print $2 "\t" $3 "\t" $4 "\t" $5 }' <<< "$trim_res")"
        r="$(awk -F'\t' -v s="$site" '$1 == s { print $2 "\t" $3 }' <<< "$raw_res")"
        b="$(awk -F'\t' -v s="$site" '$1 == s { print $2 }' "$OUT/.bam_tmp")"; b="${b:-NA}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run" "$strain" "$library" "$site" "$b" "$t" "$r" >> "$OUT/counts.tsv"
    done < "$OUT/sites.tsv"
done
rm -f "$OUT/.bam_tmp"

# ── summary to the log: pooled by strain, library, site ─────────────────────────
log_info "pooled over replicates (bam, trimmed, trimmed with A at +1, trimmed mm5, raw, raw with an all-G prefix):"
awk -F'\t' 'NR > 1 { k = $2 "\t" $3 "\t" $4; b[k] += $5; t[k] += $6; a[k] += $7; m[k] += $8; r[k] += $10; g[k] += $11 }
            END { for (k in b) print k "\t" b[k] "\t" t[k] "\t" a[k] "\t" m[k] "\t" r[k] "\t" g[k] }' "$OUT/counts.tsv" \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3n \
    | (printf "strain\tlibrary\tsite\tbam\ttrimmed\ttrim_A\ttrim_mm5\traw\traw_Gpre\n"; cat) | column -t
log_info "done"
