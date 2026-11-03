#!/usr/bin/env python
"""Stage 22 step 5: do unclipped reads start where the long-clip reads start?

A long-clip read and an M read (one with no 5' soft clip) that start at the same
5' end may be the same kind of molecule aligned two ways. For a spliced clip the
M counterpart is a read that STAR split across the same intron instead of
clipping the exon-1 piece. Three tables are written to --outdir:

  lc_m_profile.tsv   class, offset, kind, lc_reads, lc_with, pct
        For each origin class and each offset d in -W..W along the read strand
        from the long-clip read's inferred 5' end (tss5), the share of long-clip
        reads with at least one M read starting at tss5 + d. kind "any_M" counts
        every M read; kind "spliced_M_same_intron" counts only M reads whose
        5'-most intron is the long-clip read's intron (splice classes only).
        Offsets other than 0 are the local background.

  lc_m_exon1.tsv     intron_id, exon1_len, kind, n_reads
        For every intron that a splice-class long-clip read uses, the reads that
        start in exon 1 and reach exon 2 through it, by the length of their
        exon-1 piece. kind is "clipped_extra5" or "clipped_no_extra5" (a
        long-clip read, the piece being its templated clip, with or without
        untemplated 5' bases), "split_M" (an M read split across the intron) or
        "split_short_clip" (a read with a 5' clip of 1..max_clip bases, e.g.
        1Sg, split across the intron). The pairs to compare are clipped_extra5
        with split_short_clip, and clipped_no_extra5 with split_M.

  lc_m_example.tsv   read, view, start, end, part, chr, strand, tss5, intron_id
        The drawable blocks of a few long-clip and split M reads at the splice
        position with the most reads of both kinds, for the figure's panel A.
        Only long-clip reads without a clip mismatch are drawn, as the typical
        case (97 % of splice_annotated reads).
        view is "STAR alignment", "realigned" or "split M read"; part is
        soft_clip, untemplated, exon1_piece, aligned (exon 2 side) or intron.

Usage (from project root):
    python scripts/14_longclip/05_m_counterparts.py --samples config/samples.tsv \\
        --reads results/longclip/longclip_reads.tsv.gz \\
        --realigned results/longclip/bam/longclip_realigned.bam \\
        --max-clip 4 --outdir results/longclip [--window 20] [--force]
"""
import argparse
import csv
import gzip
import os
import re
from collections import Counter, defaultdict

import pysam

S, M, I, D, N, EQ, X = 4, 0, 1, 2, 3, 7, 8
REF_OPS = (M, D, N, EQ, X)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samples", required=True)
    p.add_argument("--reads", required=True)
    p.add_argument("--realigned", required=True)
    p.add_argument("--max-clip", type=int, required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--window", type=int, default=20)
    p.add_argument("--n-example", type=int, default=6)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def sample_bams(path):
    out, header = [], None
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            continue
        row = dict(zip(header, f))
        if row.get("bam") and os.path.exists(row["bam"]):
            out.append(row["bam"])
    return out


def ref_blocks(a):
    """Aligned blocks and introns of a record, 1-based inclusive, genome order."""
    r, blocks, introns = a.reference_start, [], []
    for op, n in a.cigartuples:
        if op in (M, EQ, X):
            if blocks and blocks[-1][1] == r:      # merge across =/X boundaries
                blocks[-1] = (blocks[-1][0], r + n)
            else:
                blocks.append((r + 1, r + n))
        elif op == D and blocks:
            blocks[-1] = (blocks[-1][0], r + n)
        elif op == N:
            introns.append((r + 1, r + n))
        if op in REF_OPS:
            r += n
    return blocks, introns


def main():
    a = parse_args()
    prof_f = os.path.join(a.outdir, "lc_m_profile.tsv")
    ex1_f = os.path.join(a.outdir, "lc_m_exon1.tsv")
    exm_f = os.path.join(a.outdir, "lc_m_example.tsv")
    if not a.force and all(os.path.exists(f) for f in (prof_f, ex1_f, exm_f)):
        print(f"Output exists, skipping (use --force to rerun): {prof_f}")
        return

    rows = list(csv.DictReader(gzip.open(a.reads, "rt"), delimiter="\t"))
    splice_introns = {r["intron_id"] for r in rows
                      if r["class"].startswith("splice") and r["intron_id"]}
    print(f"Long-clip rows: {len(rows)}; introns used by splice classes: {len(splice_introns)}")

    m_any = Counter()                  # (chr, strand, pos5) -> M reads
    m_intron = Counter()               # (chr, strand, pos5, intron_id) -> spliced M reads
    exon1 = Counter()                  # (intron_id, exon1_len, kind) -> reads
    split_examples = defaultdict(list)  # (chr, strand, pos5, intron_id) -> blocks
    n_m = 0
    for bam in sample_bams(a.samples):
        with pysam.AlignmentFile(bam) as fh:
            for r in fh:
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                ct = r.cigartuples
                rev = r.is_reverse
                clip5 = (ct[-1][1] if ct[-1][0] == S else 0) if rev else \
                        (ct[0][1] if ct[0][0] == S else 0)
                if clip5 > a.max_clip:
                    continue
                st = "-" if rev else "+"
                p5 = r.reference_end if rev else r.reference_start + 1
                if clip5 == 0:
                    n_m += 1
                    m_any[(r.reference_name, st, p5)] += 1
                if not any(op == N for op, _ in ct):
                    continue
                blocks, introns = ref_blocks(r)
                i5 = introns[-1] if rev else introns[0]
                iid = f"{r.reference_name}:{i5[0]}-{i5[1]}:{st}"
                if iid not in splice_introns:
                    continue
                b5 = blocks[-1] if rev else blocks[0]
                kind = "split_M" if clip5 == 0 else "split_short_clip"
                exon1[(iid, b5[1] - b5[0] + 1, kind)] += 1
                if clip5 == 0:
                    key = (r.reference_name, st, p5, iid)
                    m_intron[key] += 1
                    if len(split_examples[key]) < a.n_example:
                        split_examples[key].append((r.query_name, blocks, introns))
    print(f"M reads: {n_m}")

    # ── profile around each long-clip 5' end ──────────────────────────────────
    acc = defaultdict(lambda: [0, 0])
    for r in rows:
        cls, n = r["class"], int(r["n_reads"])
        sign = 1 if r["strand"] == "+" else -1
        t5 = int(r["tss5"])
        for d in range(-a.window, a.window + 1):
            k = (r["chr"], r["strand"], t5 + sign * d)
            acc[(cls, d, "any_M")][0] += n
            acc[(cls, d, "any_M")][1] += n * (m_any[k] > 0)
            if cls.startswith("splice"):
                acc[(cls, d, "spliced_M_same_intron")][0] += n
                acc[(cls, d, "spliced_M_same_intron")][1] += n * (m_intron[k + (r["intron_id"],)] > 0)
    with open(prof_f, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["class", "offset", "kind", "lc_reads", "lc_with", "pct"])
        for (cls, d, kind), (tot, wi) in sorted(acc.items()):
            w.writerow([cls, d, kind, tot, wi, round(100 * wi / tot, 3)])
    print(f"Written: {prof_f}")

    # ── exon-1 piece length: clipped vs split ─────────────────────────────────
    for r in rows:
        if r["class"].startswith("splice") and r["intron_id"]:
            L = int(r["tmpl_end"]) - int(r["tmpl_start"]) + 1
            nu = int(r["n_untemplated"] or 0)
            exon1[(r["intron_id"], L, "clipped_extra5" if nu else "clipped_no_extra5")] += int(r["n_reads"])
    with open(ex1_f, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["intron_id", "exon1_len", "kind", "n_reads"])
        for (iid, L, kind), n in sorted(exon1.items()):
            w.writerow([iid, L, kind, n])
    print(f"Written: {ex1_f}")

    # ── example locus ─────────────────────────────────────────────────────────
    # The splice_annotated 5' end with the most reads of both kinds. Each drawn
    # long-clip read appears twice: as STAR aligned it (soft clip beside the
    # aligned 5' end) and as step 4 realigned it (exon-1 piece, intron, exon 2).
    lc_at = Counter()
    for r in rows:
        if r["class"] == "splice_annotated":
            lc_at[(r["chr"], r["strand"], int(r["tss5"]), r["intron_id"])] += int(r["n_reads"])
    best = max(lc_at, key=lambda k: min(lc_at[k], m_intron[k]))
    chrom, st, t5, iid = best
    i_s, i_e = (int(x) for x in iid.split(":")[1].split("-"))
    print(f"Example: {iid}, 5' end {chrom}:{t5}:{st}, "
          f"{lc_at[best]} long-clip and {m_intron[best]} split M reads")
    out = []

    def exon_parts(name, label, blocks, introns, nu):
        for b in blocks:
            on_exon1 = b[0] > i_e if st == "-" else b[1] < i_s
            out.append((name, label, b[0], b[1], "exon1_piece" if on_exon1 else "aligned"))
        for i in introns:
            out.append((name, label, i[0], i[1], "intron"))
        if nu:
            u = (t5 + 1, t5 + nu) if st == "-" else (t5 - nu, t5 - 1)
            out.append((name, label, u[0], u[1], "untemplated"))

    picked, seen_oc = [], set()
    with pysam.AlignmentFile(a.realigned) as fh:
        for r in fh.fetch(chrom, i_s - 200, i_e + 200):
            if r.get_tag("lc") != "splice_annotated" or r.get_tag("cm") \
                    or ("-" if r.is_reverse else "+") != st:
                continue
            blocks, introns = ref_blocks(r)
            p5 = blocks[-1][1] if r.is_reverse else blocks[0][0]
            if p5 != t5 or not introns or (introns[-1] if r.is_reverse else introns[0]) != (i_s, i_e):
                continue
            key = (r.get_tag("oc"), r.get_tag("op"))
            if key in seen_oc:
                continue
            seen_oc.add(key)
            picked.append(r)
    picked = picked[:a.n_example // 2 or 1]
    for n, r in enumerate(picked, 1):
        name = f"long-clip {n}"
        oc, op_ = r.get_tag("oc"), r.get_tag("op")
        # STAR's alignment, rebuilt from the original CIGAR and POS
        pos, blocks = op_, []
        for ln, o in re.findall(r"(\d+)([MIDNSHP=X])", oc):
            ln = int(ln)
            if o in "M=X":
                blocks.append((pos, pos + ln - 1)); pos += ln
            elif o in "DN":
                pos += ln
        clen = int(re.findall(r"(\d+)S", oc)[-1 if st == "-" else 0])
        lo, hi = blocks[0][0], blocks[-1][1]
        sc = (hi + 1, hi + clen) if st == "-" else (lo - clen, lo - 1)
        for b in blocks:
            out.append((name, "STAR alignment", b[0], b[1], "aligned"))
        out.append((name, "STAR alignment", sc[0], sc[1], "soft_clip"))
        rb, ri = ref_blocks(r)
        exon_parts(name, "realigned", rb, ri, r.get_tag("nu"))
    for k, (q, blocks, introns) in enumerate(split_examples[best][:a.n_example // 2 or 1], 1):
        exon_parts(f"split M {k}", "split M read", blocks, introns, 0)
    with open(exm_f, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["read", "view", "start", "end", "part", "chr", "strand", "tss5", "intron_id"])
        for row in out:
            w.writerow(list(row) + [chrom, st, t5, iid])
    print(f"Written: {exm_f}")


if __name__ == "__main__":
    main()
