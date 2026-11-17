#!/usr/bin/env python
"""NAD branch step 9, part 1: per-position fine 5' clip types from CAGE BAMs.

For every uniquely mapped primary alignment (MAPQ 255) the 5' end position is
taken strand-aware (POS for +, reference end for -), and the 5' soft clip is
written 5'->3' in RNA sense (reverse-complemented for - strand reads), so the
last letter of the label is the base that abuts the genome, as in pyselectal.
Clips of 1-3 nt keep their full sequence (1Sa, 2Sta, 3Sgta ...), longer clips
are pooled as "4+S", and reads without a 5' clip are "M". Unlike
results/tss/typed_ctss.tsv.gz, nothing is collapsed into "other", so types such
as 2Sta and 2Sca stay visible.

Usage (from project root, in an env with pysam):
  python nad/scripts/09_chrM_clip_spectrum/fine_clip_types.py OUT.tsv.gz BAM [BAM ...]
"""
import collections
import gzip
import sys

import pysam

COMP = str.maketrans("ACGTN", "TGCAN")


def main(out_path, bams):
    counts = collections.Counter()
    for bam in bams:
        with pysam.AlignmentFile(bam) as fh:
            for r in fh:
                if (r.is_unmapped or r.is_secondary or r.is_supplementary
                        or r.mapping_quality != 255):
                    continue
                cig, seq = r.cigartuples, r.query_sequence
                if r.is_reverse:
                    op, n = cig[-1]
                    pos, strand = r.reference_end, "-"
                    clip = seq[-n:].translate(COMP)[::-1] if op == 4 else ""
                else:
                    op, n = cig[0]
                    pos, strand = r.reference_start + 1, "+"
                    clip = seq[:n] if op == 4 else ""
                if not clip:
                    t = "M"
                elif len(clip) <= 3:
                    t = f"{len(clip)}S{clip.lower()}"
                else:
                    t = "4+S"
                counts[(r.reference_name, pos, strand, t)] += 1
    with gzip.open(out_path, "wt") as f:
        f.write("chr\tpos\tstrand\ttype\tn\n")
        for (c, p, st, t), n in sorted(counts.items()):
            f.write(f"{c}\t{p}\t{st}\t{t}\t{n}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
