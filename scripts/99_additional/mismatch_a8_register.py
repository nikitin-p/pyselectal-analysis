#!/usr/bin/env python3
"""Tally the upstream base composition of reads by their first mismatch position.

The question is whether a read of the M class that mismatches the genome at its
third aligned base really starts at its first base. The alternative is that the
mismatching base is an untemplated addition sitting at -1 of the true initiator,
so that the true +1 is the fourth aligned base and the first two bases matched
the genome by chance. The yeast A at -8 serves as the ruler: it sits 8 bases
upstream of a real TSS, so it marks the register in which a group of reads
starts.

Every read is assigned to a group by the position of its first mismatch, counted
from the 5' end of the alignment in the orientation of the read:
    none      no mismatch within the first --down aligned bases
    3 .. 8    the first mismatch is at that position
    3G 3A 3Y  the same reads as group 3, split by the mismatching read base
    late      the first mismatch lies beyond position --kmax
Reads with an insertion, a deletion or a splice are skipped, because their read
positions do not map one to one onto genomic offsets.

Three tables are written to --outdir:
    a8_register_profile.tsv  type, group, weight, coord, base, n
    a8_register_dinuc.tsv    type, group, ref_k, read_k, ref_next, n
    a8_register_reads.tsv    type, group, n_reads, n_sites
The coord column has no zero: +1 is the first aligned base and -1 is the base
before it. The weight column is "reads" when every read counts and "sites" when
every distinct 5' position counts once.

Usage (from project root):
    python scripts/99_additional/mismatch_a8_register.py --types M,1Sg \\
        [--samples config/samples.tsv] [--bam_dir results/bam_selected] \\
        [--fasta data/sacCer3.fa] [--outdir results/mismatch] \\
        [--up 25] [--down 12] [--kmax 8] [--force]
"""
import argparse
import collections
import os
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
OUTPUTS = ("a8_register_profile.tsv", "a8_register_dinuc.tsv",
           "a8_register_reads.tsv")


def read_sample_ids(path):
    ids = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index("sample_id")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ids.append(line.rstrip("\n").split("\t")[col])
    return ids


def window(chrom, pos5, reverse, up, down):
    """Return the genomic window around a 5' position, in read orientation."""
    if reverse:
        lo, hi = pos5 + 1 - down, pos5 + 1 + up
        if lo < 0 or hi > len(chrom):
            return None
        return chrom[lo:hi].translate(COMPLEMENT)[::-1]
    lo, hi = pos5 - up, pos5 + down
    if lo < 0 or hi > len(chrom):
        return None
    return chrom[lo:hi]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--types", default="M",
                    help="comma-separated end types; one BAM per sample each")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--bam_dir", default="results/bam_selected")
    ap.add_argument("--fasta", default="data/sacCer3.fa")
    ap.add_argument("--outdir", default="results/mismatch")
    ap.add_argument("--up", type=int, default=25)
    ap.add_argument("--down", type=int, default=12)
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    up, down = args.up, args.down

    out_paths = [os.path.join(args.outdir, f) for f in OUTPUTS]
    if all(os.path.exists(p) for p in out_paths) and not args.force:
        print(f"Outputs exist in {args.outdir}; skipping (use --force).",
              file=sys.stderr)
        return
    os.makedirs(args.outdir, exist_ok=True)

    with pysam.FastaFile(args.fasta) as fa:
        ref = {c: fa.fetch(c).upper() for c in fa.references}
    samples = read_sample_ids(args.samples)

    sites = collections.defaultdict(collections.Counter)  # (type, group) -> site -> reads
    dinuc = collections.Counter()   # (type, group, ref_k, read_k, ref_next) -> reads
    for t in args.types.split(","):
        for s in samples:
            bam_path = os.path.join(args.bam_dir, f"{s}_{t}.bam")
            if not os.path.exists(bam_path):
                sys.exit(f"BAM not found: {bam_path}")
            with pysam.AlignmentFile(bam_path) as bam:
                for r in bam:
                    if any(op in (1, 2, 3) for op, _ in r.cigartuples):
                        continue
                    seq = r.query_sequence
                    if r.is_reverse:
                        pos5 = r.reference_end - 1
                        rd = seq[:r.query_alignment_end].translate(COMPLEMENT)[::-1]
                    else:
                        pos5 = r.reference_start
                        rd = seq[r.query_alignment_start:]
                    rd = rd[:down]
                    win = window(ref[r.reference_name], pos5, r.is_reverse, up, down)
                    if win is None or len(rd) < down:
                        continue
                    g = win[up:]
                    first = next((i + 1 for i in range(down) if rd[i] != g[i]), 0)
                    if first == 0:
                        groups = ["none"]
                    elif first > args.kmax:
                        groups = ["late"]
                    else:
                        groups = [str(first)]
                        dinuc[(t, first, g[first - 1], rd[first - 1], g[first])] += 1
                        if first == 3:
                            base = rd[2] if rd[2] in "GA" else "Y"
                            groups.append("3" + base)
                    site = (r.reference_name, r.is_reverse, pos5)
                    for grp in groups:
                        sites[(t, grp)][site] += 1
            print(f"{t}\t{s}", file=sys.stderr)

    coords = [i - up if i < up else i - up + 1 for i in range(up + down)]
    with open(out_paths[0], "w") as prof, open(out_paths[2], "w") as nfh:
        prof.write("type\tgroup\tweight\tcoord\tbase\tn\n")
        nfh.write("type\tgroup\tn_reads\tn_sites\n")
        for (t, grp), counter in sorted(sites.items()):
            by_reads = [collections.Counter() for _ in coords]
            by_sites = [collections.Counter() for _ in coords]
            for (chrom, reverse, pos5), n in counter.items():
                for i, b in enumerate(window(ref[chrom], pos5, reverse, up, down)):
                    by_reads[i][b] += n
                    by_sites[i][b] += 1
            for weight, tally in (("reads", by_reads), ("sites", by_sites)):
                for i, c in enumerate(coords):
                    for b in sorted(tally[i]):
                        prof.write(f"{t}\t{grp}\t{weight}\t{c}\t{b}\t{tally[i][b]}\n")
            nfh.write(f"{t}\t{grp}\t{sum(counter.values())}\t{len(counter)}\n")
    with open(out_paths[1], "w") as fh:
        fh.write("type\tgroup\tref_k\tread_k\tref_next\tn\n")
        for k in sorted(dinuc):
            fh.write("\t".join(map(str, k + (dinuc[k],))) + "\n")
    for p in out_paths:
        print(f"Wrote {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
