"""NAD branch step 4: count strand-aware 5' ends of NAD captureSeq alignments.

Reads are sense to the RNA, so the 5' end is the leftmost aligned base on the
+ strand and the rightmost on the - strand. Positions are 1-based, matching
results/tss/typed_ctss.tsv.gz (Stage 8, GenomicAlignments start/end).

Read names carry  <id>_g<Glen>_<UMI>  from step 2. Per position the script
reports:
  reads      every unique alignment;
  umis       UMI-collapsed count: distinct UMIs among reads with a known UMI,
             plus every read whose UMI is unknown (NNNNNN), which cannot be
             collapsed;
  mm5        reads whose first aligned base mismatches the genome, i.e. a
             non-templated 5' base (alignment is Extend5pOfRead1, no clip);
  g0         reads with no leading G stripped in step 2 (useful because a
             genomic G at +1 is ambiguous with the TdT tail).
The genomic base at +1 (strand-aware) is added as base_p1.

Usage:
  python count_5p.py --bam RUN.bam --fasta sacCer3.fa --out RUN.ctss.tsv.gz
"""
import argparse
import gzip
from collections import defaultdict

import pysam

COMP = str.maketrans("ACGTN", "TGCAN")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mapq", type=int, default=255)
    a = ap.parse_args()

    fa = pysam.FastaFile(a.fasta)
    reads = defaultdict(int)
    umiset = defaultdict(set)
    nounk = defaultdict(int)
    mm5 = defaultdict(int)
    g0 = defaultdict(int)

    with pysam.AlignmentFile(a.bam) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.mapq:
                continue
            parts = r.query_name.rsplit("_", 2)
            umi, glen = parts[-1], parts[-2]
            seq = r.query_sequence
            chrom = r.reference_name
            if r.is_reverse:
                p0 = r.reference_end - 1          # 0-based
                key = (chrom, p0 + 1, "-")
                read_base = seq[-1]
            else:
                p0 = r.reference_start
                key = (chrom, p0 + 1, "+")
                read_base = seq[0]
            reads[key] += 1
            if umi == "NNNNNN" or "N" in umi:
                nounk[key] += 1
            else:
                umiset[key].add(umi)
            if read_base != fa.fetch(chrom, p0, p0 + 1).upper():
                mm5[key] += 1
            if glen == "g0":
                g0[key] += 1

    with gzip.open(a.out, "wt") as out:
        out.write("chr\tpos\tstrand\tbase_p1\treads\tumis\tmm5\tg0\n")
        for key in sorted(reads, key=lambda k: (k[0], k[1], k[2])):
            chrom, pos, strand = key
            b = fa.fetch(chrom, pos - 1, pos).upper()
            if strand == "-":
                b = b.translate(COMP)
            out.write(f"{chrom}\t{pos}\t{strand}\t{b}\t{reads[key]}\t"
                      f"{len(umiset[key]) + nounk[key]}\t{mm5[key]}\t{g0[key]}\n")


if __name__ == "__main__":
    main()
