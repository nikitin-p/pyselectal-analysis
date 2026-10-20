#!/usr/bin/env python3
"""Profile the mismatches of per-type BAMs against the genome, by read position.

The question is whether the reads of the M class, which carry no 5' soft clip,
match the genome throughout or only at their 5' end. Every aligned base of
every read is compared with the reference, and each mismatch is recorded at
its position counted from the 5' end of the alignment, in the orientation of
the read, so that position 1 is the first aligned base on either strand. A
second type given in --types (normally the canonical one) is profiled in the
same way and serves as the baseline in mismatch_plot.R.

Under STAR --alignEndsType Local a mismatch at position 1 or 2 is soft-clipped
rather than aligned, because aligning it would not raise the score. Those two
positions are therefore expected to be free of mismatches by construction, and
position 3 is the first one at which a mismatch can survive.

The BAMs carry NM but not MD, so the reference is read from the FASTA. The
comparison walks the CIGAR by hand: M blocks are compared, I and S advance the
read only, D and N advance the reference only.

Five tables are written to --outdir:
    position_substitutions.tsv   type, sample_id, pos, ref, read, qual, n
    position_coverage.tsv        type, sample_id, pos, n
    reads_by_mismatch_count.tsv  type, sample_id, n_mismatch, has_indel, n
    prefix3.tsv                  type, sample_id, prefix, n_reads, n_pos3_mismatch
    prefix3_ref.tsv              type, sample_id, prefix, ref, n
The qual column is "hi" for a base quality of at least --min_qual and "lo"
otherwise. The prefix is the first three aligned bases of the read, 5'->3'.
prefix3_ref.tsv splits the position-3 mismatches of each prefix by the genomic
base they replace, in the orientation of the read.

Usage (from project root):
    python scripts/99_additional/mismatch_profile.py --types M,1Sg \\
        [--samples config/samples.tsv] [--bam_dir results/bam_selected] \\
        [--fasta data/sacCer3.fa] [--outdir results/mismatch] \\
        [--min_qual 30] [--force]
"""
import argparse
import collections
import os
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
OUTPUTS = ("position_substitutions.tsv", "position_coverage.tsv",
           "reads_by_mismatch_count.tsv", "prefix3.tsv", "prefix3_ref.tsv")


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


def profile(bam_path, ref, min_qual):
    """Return the five counters for one BAM."""
    subs = collections.Counter()      # (pos, ref base, read base, qual class)
    alen_n = collections.Counter()    # aligned query length -> reads
    reads = collections.Counter()     # (mismatches, has indel) -> reads
    prefix = collections.Counter()    # 3-mer -> reads
    prefix_mm = collections.Counter() # 3-mer -> reads mismatching at position 3
    prefix_ref = collections.Counter() # (3-mer, genomic base at position 3) -> reads
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam:
            seq, qual = r.query_sequence, r.query_qualities
            chrom = ref[r.reference_name]
            qs, qe = r.query_alignment_start, r.query_alignment_end
            rp, qp = r.reference_start, 0
            mismatches, has_indel = [], False
            for op, length in r.cigartuples:
                if op in (0, 7, 8):
                    a, b = seq[qp:qp + length], chrom[rp:rp + length]
                    if a != b:
                        for i in range(length):
                            if a[i] != b[i]:
                                mismatches.append((qp + i, b[i], a[i]))
                    qp += length
                    rp += length
                elif op == 4:
                    qp += length
                elif op == 1:
                    has_indel = True
                    qp += length
                elif op == 2:
                    has_indel = True
                    rp += length
                elif op == 3:
                    rp += length
            ref3 = None
            for qi, rb, qb in mismatches:
                if r.is_reverse:
                    pos = qe - qi
                    rb, qb = rb.translate(COMPLEMENT), qb.translate(COMPLEMENT)
                else:
                    pos = qi - qs + 1
                if pos == 3:
                    ref3 = rb
                subs[(pos, rb, qb, "hi" if qual[qi] >= min_qual else "lo")] += 1
            alen_n[qe - qs] += 1
            reads[(len(mismatches), int(has_indel))] += 1
            if qe - qs >= 3:
                if r.is_reverse:
                    p3 = seq[qe - 3:qe].translate(COMPLEMENT)[::-1]
                else:
                    p3 = seq[qs:qs + 3]
                prefix[p3] += 1
                if ref3 is not None:
                    prefix_mm[p3] += 1
                    prefix_ref[(p3, ref3)] += 1
    # A read of aligned length L covers positions 1..L.
    coverage, running = {}, 0
    for pos in range(max(alen_n, default=0), 0, -1):
        running += alen_n.get(pos, 0)
        coverage[pos] = running
    return subs, coverage, reads, prefix, prefix_mm, prefix_ref


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--types", default="M",
                    help="comma-separated end types; one BAM per sample each")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--bam_dir", default="results/bam_selected")
    ap.add_argument("--fasta", default="data/sacCer3.fa")
    ap.add_argument("--outdir", default="results/mismatch")
    ap.add_argument("--min_qual", type=int, default=30)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_paths = [os.path.join(args.outdir, f) for f in OUTPUTS]
    if all(os.path.exists(p) for p in out_paths) and not args.force:
        print(f"Outputs exist in {args.outdir}; skipping (use --force).",
              file=sys.stderr)
        return
    os.makedirs(args.outdir, exist_ok=True)

    with pysam.FastaFile(args.fasta) as fa:
        ref = {c: fa.fetch(c).upper() for c in fa.references}
    samples = read_sample_ids(args.samples)
    types = args.types.split(",")

    rows = {f: [] for f in OUTPUTS}
    for t in types:
        for s in samples:
            bam_path = os.path.join(args.bam_dir, f"{s}_{t}.bam")
            if not os.path.exists(bam_path):
                sys.exit(f"BAM not found: {bam_path}")
            subs, coverage, reads, prefix, prefix_mm, prefix_ref = profile(
                bam_path, ref, args.min_qual)
            for k in sorted(subs):
                rows[OUTPUTS[0]].append((t, s) + k + (subs[k],))
            for pos in sorted(coverage):
                rows[OUTPUTS[1]].append((t, s, pos, coverage[pos]))
            for k in sorted(reads):
                rows[OUTPUTS[2]].append((t, s) + k + (reads[k],))
            for k in sorted(prefix):
                rows[OUTPUTS[3]].append((t, s, k, prefix[k], prefix_mm[k]))
            for k in sorted(prefix_ref):
                rows[OUTPUTS[4]].append((t, s) + k + (prefix_ref[k],))
            print(f"{t}\t{s}\t{sum(reads.values())} reads", file=sys.stderr)

    headers = {
        OUTPUTS[0]: ("type", "sample_id", "pos", "ref", "read", "qual", "n"),
        OUTPUTS[1]: ("type", "sample_id", "pos", "n"),
        OUTPUTS[2]: ("type", "sample_id", "n_mismatch", "has_indel", "n"),
        OUTPUTS[3]: ("type", "sample_id", "prefix", "n_reads", "n_pos3_mismatch"),
        OUTPUTS[4]: ("type", "sample_id", "prefix", "ref", "n"),
    }
    for f, path in zip(OUTPUTS, out_paths):
        with open(path, "w") as fh:
            fh.write("\t".join(headers[f]) + "\n")
            for row in rows[f]:
                fh.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
