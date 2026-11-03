#!/usr/bin/env python
"""Stage 22 step 4: every long-clip read as a match/mismatch BAM record.

Each long-clip read of the subsampled BAMs is joined to its Stage 22 row on
<chr, strand, aln_start, aln_end, cigar, clip> and realigned so that its whole
sequence, 5' to 3', is compared with the genome:

  * a placed clip (splice_*, templated_*) is aligned at tmpl_start..tmpl_end,
    joined to the aligned part by an N gap (D for templated_adjacent), with its
    n_untemplated leading bases kept soft-clipped;
  * an unplaced clip is aligned onto the clen bases directly upstream of the
    aligned 5' end;
  * a 3' soft clip is aligned onto the bases directly downstream;
  * a clip that would run off the chromosome stays soft-clipped.

Every aligned base becomes '=' or 'X' in the CIGAR. BAM stores bases in four
bits and loses case, so the case-encoded read (UPPER match, lower mismatch,
insertion, untemplated or unplaceable clip base), written 5'->3' in the read's
orientation, goes into the ms:Z tag and the per-read table.

Outputs in --outdir:
    longclip_realigned.bam(.bai)  one record per read; tags ms YC lc nu cm am oc op RG
    longclip_segments.bam(.bai)   one record per maximal match / mismatch /
                                  untemplated run; BX:Z = read name, YC = colour
    longclip_case.tsv.gz          one row per read

Usage (from project root):
    python scripts/14_longclip/04_case_bams.py --samples config/samples.tsv \\
        --reads results/longclip/longclip_reads.tsv.gz --fasta data/sacCer3.fa \\
        --helpers scripts/14_longclip/lc_helpers.R --max-clip 4 \\
        --outdir results/longclip/bam [--force]
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

import pysam

S, M, I, D, N = 4, 0, 1, 2, 3
EQ, X = 7, 8
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

MATCH_RGB, MISMATCH_RGB, UNTEMPL_RGB = "0,150,0", "220,38,38", "128,128,128"


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samples", required=True)
    p.add_argument("--reads", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--helpers", required=True,
                   help="lc_helpers.R, read for the origin-class palette")
    p.add_argument("--max-clip", type=int, required=True,
                   help="dom_max_clip; a clip longer than this is a long clip")
    p.add_argument("--outdir", required=True)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def class_colours(helpers):
    """LC_CLASS_COLORS from lc_helpers.R, as 'r,g,b' strings."""
    txt = open(helpers).read()
    block = re.search(r"LC_CLASS_COLORS\s*<-\s*c\((.*?)\)", txt, re.S).group(1)
    out = {}
    for name, hx in re.findall(r"(\w+)\s*=\s*\"#([0-9a-fA-F]{6})\"", block):
        out[name] = ",".join(str(int(hx[i:i + 2], 16)) for i in (0, 2, 4))
    return out


def read_samples(path):
    out = []
    with open(path) as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            row = dict(zip(header, f))
            if row.get("bam") and os.path.exists(row["bam"]):
                out.append((row["sample_id"], row["bam"]))
            else:
                print(f"BAM not found, skipping: {row.get('bam')}", file=sys.stderr)
    return out


def ival(x):
    return int(x) if x not in ("", "NA") else None


def load_reads(path):
    table = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["chr"], r["strand"], int(r["aln_start"]), int(r["aln_end"]),
                   r["cigar"], r["clip"])
            table[key] = {
                "read_key": r["read_key"], "class": r["class"],
                "nu": ival(r["n_untemplated"]) or 0,
                "tmpl_start": ival(r["tmpl_start"]), "tmpl_end": ival(r["tmpl_end"]),
                "gap": ival(r["gap"]), "n_reads": int(r["n_reads"]),
                "intron": tuple(int(x) for x in
                                r["intron_id"].split(":")[1].split("-")) if r["intron_id"] else None,
                "identity": float(r["identity"]) if r["identity"] not in ("", "NA") else None,
            }
    return table


def trim_ref(core, k, from_start):
    """Drop k reference bases from one end of the core CIGAR, which must be pure
    M there; return the rest, or None when an indel or gap is in the way."""
    ops = list(core) if from_start else list(reversed(core))
    while k > 0:
        if not ops or ops[0][0] != M:
            return None
        op, n = ops.pop(0)
        if n > k:
            ops.insert(0, (op, n - k))
        k -= min(n, k)
    return ops if from_start else list(reversed(ops))


def intron_layout(row, strand, core, aln_start0, aln_end, t):
    """Place a spliced clip on its intron (the annotated one, or Stage 22's slid
    GT..AG one). STAR often extends the aligned part k bases into the intron's
    3' end, which equals the tail of exon 1; those k bases are moved back to
    exon 1 so that the N gap is the intron itself. Returns (k, trimmed core,
    intron length), or None when the read does not fit that layout."""
    if row["intron"] is None:
        return None
    i_s, i_e = row["intron"]
    if strand == "+":
        k = i_e - aln_start0                      # aligned bases inside the intron
        if k < 0 or row["tmpl_end"] + k != i_s - 1:
            return None
    else:
        k = aln_end - i_s + 1
        if k < 0 or row["tmpl_start"] - k != i_e + 1:
            return None
    rest = trim_ref(core, k, from_start=(strand == "+"))
    return None if rest is None else (k, rest, i_e - i_s + 1)


def realign(a, row, strand, ref, chrlen):
    """Return (pos0, cigar ops with per-base =/X, per-query info) for one read.

    Ops are built in genome-forward order as (op, length, region) blocks with
    region in {clip5, aln, clip3}; aligned blocks are later split into =/X.
    """
    ct = a.cigartuples
    lead = ct[0][1] if ct[0][0] == S else 0
    trail = ct[-1][1] if ct[-1][0] == S else 0
    core = ct[(1 if lead else 0):(len(ct) - 1 if trail else len(ct))]
    aln_start0, aln_end = a.reference_start, a.reference_end  # 0-based, half-open
    clip5_len = lead if strand == "+" else trail
    clip3_len = trail if strand == "+" else lead
    cls, nu = row["class"], row["nu"]
    placed = cls != "unplaced"
    gap_op = D if cls == "templated_adjacent" else N
    if placed:
        # splice_annotated rows carry no gap, so it is derived from the
        # coordinates for every class and checked against the table where set.
        gap = (aln_start0 - row["tmpl_end"]) if strand == "+" else \
              (row["tmpl_start"] - aln_end - 1)
        assert gap >= 0 and (row["gap"] is None or row["gap"] == gap), row
    lay = intron_layout(row, strand, core, aln_start0, aln_end, clip5_len - nu) \
        if cls.startswith("splice") else None
    core_regions = [(op, n, "aln") for op, n in core]

    left, right = [], []
    pos0 = aln_start0
    # 5' clip
    if strand == "+":
        if placed:
            t = clip5_len - nu
            assert row["tmpl_end"] - row["tmpl_start"] + 1 == t, row
            if lay:
                k, rest, ilen = lay
                left = ([(S, nu, "clip5u")] if nu else []) + [(M, t, "clip5")] + \
                       ([(M, k, "aln")] if k else []) + [(N, ilen, None)]
                core_regions = [(op, n, "aln") for op, n in rest]
            else:
                left = ([(S, nu, "clip5u")] if nu else []) + [(M, t, "clip5")] + \
                       ([(gap_op, gap, None)] if gap else [])
            pos0 = row["tmpl_start"] - 1
        elif aln_start0 - clip5_len >= 0:
            left, pos0 = [(M, clip5_len, "clip5")], aln_start0 - clip5_len
        else:
            left = [(S, clip5_len, "clip5s")]
        if clip3_len:
            right = [(M, clip3_len, "clip3")] if aln_end + clip3_len <= chrlen \
                else [(S, clip3_len, "clip3s")]
    else:
        if placed:
            t = clip5_len - nu
            assert row["tmpl_end"] - row["tmpl_start"] + 1 == t, row
            if lay:
                k, rest, ilen = lay
                right = [(N, ilen, None)] + ([(M, k, "aln")] if k else []) + \
                        [(M, t, "clip5")] + ([(S, nu, "clip5u")] if nu else [])
                core_regions = [(op, n, "aln") for op, n in rest]
            else:
                right = ([(gap_op, gap, None)] if gap else []) + \
                        [(M, t, "clip5")] + ([(S, nu, "clip5u")] if nu else [])
        elif aln_end + clip5_len <= chrlen:
            right = [(M, clip5_len, "clip5")]
        else:
            right = [(S, clip5_len, "clip5s")]
        if clip3_len:
            if aln_start0 - clip3_len >= 0:
                left, pos0 = [(M, clip3_len, "clip3")], aln_start0 - clip3_len
            else:
                left = [(S, clip3_len, "clip3s")]
    blocks = left + core_regions + right

    seq, qual = a.query_sequence, a.query_qualities
    q, r = 0, pos0
    ops = []                 # expanded, one entry per base or per D/N block
    per_q = [None] * len(seq)  # (state, region, ref pos0 or None)
    for op, n, region in blocks:
        if op in (M, EQ, X):
            for _ in range(n):
                st = EQ if seq[q].upper() == ref[r] and ref[r] != "N" else X
                ops.append((st, 1))
                per_q[q] = (st, region, r)
                q += 1; r += 1
        elif op in (S, I):
            for _ in range(n):
                ops.append((op, 1))
                per_q[q] = (op, region, None)
                q += 1
        elif op in (D, N):
            ops.append((op, n)); r += n
        else:
            raise ValueError(f"unexpected CIGAR op {op} in {a.query_name}")
    assert q == len(seq), a.query_name
    rle = []
    for op, n in ops:
        if rle and rle[-1][0] == op:
            rle[-1] = (op, rle[-1][1] + n)
        else:
            rle.append((op, n))
    return pos0, rle, per_q, lay is not None


def case_string(seq, per_q, strand):
    chars = [b.upper() if st == EQ else b.lower() for b, (st, _, _) in zip(seq, per_q)]
    s = "".join(chars)
    return s if strand == "+" else revcomp(s)


def longest_run(states):
    best = cur = 0
    for st in states:
        cur = cur + 1 if st == EQ else 0
        best = max(best, cur)
    return best


def segments(a, per_q, strand, row):
    """Maximal runs of match / mismatch / untemplated bases with contiguous
    reference positions, in genome order, as (state, pos0, qstart, qend)."""
    q_items = list(per_q)
    # Untemplated 5' bases have no genomic position; draw them on the bases
    # directly beyond the templated clip so that they stay visible.
    utq = [i for i, (st, reg, _) in enumerate(q_items) if reg == "clip5u"]
    if utq:
        ref_q = [i for i, (_, reg, rp) in enumerate(q_items) if reg == "clip5"]
        if strand == "+":
            anchor = q_items[min(ref_q)][2]
            for k, i in enumerate(reversed(utq)):
                q_items[i] = ("U", "clip5u", anchor - 1 - k)
        else:
            anchor = q_items[max(ref_q)][2]
            for k, i in enumerate(utq):
                q_items[i] = ("U", "clip5u", anchor + 1 + k)
    segs = []
    for i, (st, reg, rp) in enumerate(q_items):
        if rp is None or rp < 0:
            continue
        s = "U" if st == "U" else ("m" if st == EQ else "x")
        if segs and segs[-1][0] == s and segs[-1][3] == i and segs[-1][4] == rp:
            segs[-1][3] = i + 1; segs[-1][4] = rp + 1
        else:
            segs.append([s, rp, i, i + 1, rp + 1])
    return segs


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    real_bam = os.path.join(args.outdir, "longclip_realigned.bam")
    seg_bam = os.path.join(args.outdir, "longclip_segments.bam")
    tsv = os.path.join(args.outdir, "longclip_case.tsv.gz")
    if not args.force and all(os.path.exists(f) for f in (real_bam, seg_bam, tsv)):
        print(f"Output exists, skipping (use --force to rerun): {real_bam}")
        return

    colours = class_colours(args.helpers)
    table = load_reads(args.reads)
    expected = sum(v["n_reads"] for v in table.values())
    samples = read_samples(args.samples)
    fa = pysam.FastaFile(args.fasta)
    genome = {c: fa.fetch(c).upper() for c in fa.references}

    with pysam.AlignmentFile(samples[0][1]) as b:
        hdr = b.header.to_dict()
    hdr["RG"] = [{"ID": sid, "SM": sid} for sid, _ in samples]
    hdr["PG"] = hdr.get("PG", []) + [{"ID": "longclip_case_bams", "PN": "04_case_bams.py"}]
    hdr["HD"] = {"VN": "1.6", "SO": "unsorted"}
    tmp_real, tmp_seg = real_bam + ".unsorted", seg_bam + ".unsorted"
    out_r = pysam.AlignmentFile(tmp_real, "wb", header=hdr)
    out_s = pysam.AlignmentFile(tmp_seg, "wb", header=hdr)

    seen = Counter()
    per_class = Counter()
    bad_x = bad_len = bad_tmpl = bad_ident = 0
    n_seg = 0
    n_on_intron = Counter()
    with gzip.open(tsv, "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["qname", "sample_id", "read_key", "class", "chr", "strand", "pos5",
                    "clen", "nu", "cm", "am", "clip3_len", "max_run_clip", "ms"])
        for sid, path in samples:
            n_sample = 0
            with pysam.AlignmentFile(path) as bam:
                for a in bam:
                    if a.is_unmapped or a.is_secondary or a.is_supplementary:
                        continue
                    strand = "-" if a.is_reverse else "+"
                    ct = a.cigartuples
                    clen = (ct[0][1] if ct[0][0] == S else 0) if strand == "+" else \
                           (ct[-1][1] if ct[-1][0] == S else 0)
                    if clen <= args.max_clip:
                        continue
                    seq = a.query_sequence
                    clip = seq[:clen] if strand == "+" else revcomp(seq[-clen:])
                    if "N" in clip:  # Stage 22 step 1 drops ambiguous clips too
                        continue
                    key = (a.reference_name, strand, a.reference_start + 1,
                           a.reference_end, a.cigarstring, clip)
                    row = table.get(key)
                    if row is None:
                        sys.exit(f"Read not in the Stage 22 table: {a.query_name} {key}")
                    seen[key] += 1
                    n_sample += 1
                    ref = genome[a.reference_name]
                    pos0, rle, per_q, on_intron = realign(a, row, strand, ref, len(ref))
                    if row["class"].startswith("splice"):
                        n_on_intron[on_intron] += 1

                    cm = sum(1 for st, reg, _ in per_q if reg == "clip5" and st == X)
                    am = sum(1 for st, reg, _ in per_q if reg == "aln" and st == X)
                    clip3 = sum(1 for _, reg, _ in per_q if reg in ("clip3", "clip3s"))
                    clip_states = [st for st, reg, _ in per_q if reg == "clip5"]
                    if strand == "-":
                        clip_states = clip_states[::-1]
                    ms = case_string(seq, per_q, strand)
                    pos5 = a.reference_start + 1 if strand == "+" else a.reference_end

                    new = pysam.AlignedSegment(out_r.header)
                    new.query_name = a.query_name
                    new.flag = a.flag
                    new.reference_id = out_r.header.get_tid(a.reference_name)
                    new.reference_start = pos0
                    new.mapping_quality = a.mapping_quality
                    new.cigartuples = rle
                    new.query_sequence = seq
                    new.query_qualities = a.query_qualities
                    new.set_tags([("RG", sid), ("ms", ms), ("YC", colours[row["class"]]),
                                  ("lc", row["class"]), ("nu", row["nu"]),
                                  ("cm", cm), ("am", am), ("oc", a.cigarstring),
                                  ("op", a.reference_start + 1)])
                    # independent checks
                    if new.infer_query_length() != len(seq):
                        bad_len += 1
                    nx = sum(1 for qi, rp in new.get_aligned_pairs(matches_only=True)
                             if seq[qi].upper() != ref[rp] or ref[rp] == "N")
                    if nx != cm + am + sum(1 for st, reg, _ in per_q
                                           if reg == "clip3" and st == X):
                        bad_x += 1
                    if row["class"].startswith("templated") and cm:
                        bad_tmpl += 1
                    if row["class"] == "splice_annotated" and row["identity"] is not None:
                        t = len(clip_states)
                        if (t - cm) / t + 1e-9 < row["identity"]:
                            bad_ident += 1
                    out_r.write(new)
                    per_class[row["class"]] += 1

                    for k, (st, rp, qs, qe, _) in enumerate(segments(a, per_q, strand, row)):
                        sg = pysam.AlignedSegment(out_s.header)
                        sg.query_name = f"{a.query_name}.{k}"
                        sg.flag = 16 if a.is_reverse else 0
                        sg.reference_id = new.reference_id
                        sg.reference_start = rp
                        sg.mapping_quality = a.mapping_quality
                        sg.cigartuples = [(M, qe - qs)]
                        sg.query_sequence = seq[qs:qe]
                        if a.query_qualities is not None:
                            sg.query_qualities = a.query_qualities[qs:qe]
                        rgb = {"m": MATCH_RGB, "x": MISMATCH_RGB, "U": UNTEMPL_RGB}[st]
                        sg.set_tags([("RG", sid), ("BX", a.query_name), ("YC", rgb),
                                     ("sg", {"m": "match", "x": "mismatch",
                                             "U": "untemplated"}[st]),
                                     ("lc", row["class"])])
                        out_s.write(sg)
                        n_seg += 1

                    w.writerow([a.query_name, sid, row["read_key"], row["class"],
                                a.reference_name, strand, pos5, clen, row["nu"], cm, am,
                                clip3, longest_run(clip_states), ms])
            print(f"{sid}: {n_sample} long-clip reads")
    out_r.close(); out_s.close()

    total = sum(seen.values())
    print(f"Long-clip reads written: {total} (Stage 22 table: {expected})")
    print(f"Segments written: {n_seg}")
    print("Reads per class:", dict(sorted(per_class.items())))
    print(f"Spliced reads placed on their intron: {n_on_intron[True]}, "
          f"left on the template gap: {n_on_intron[False]}")
    off = [k for k, v in table.items() if seen[k] != v["n_reads"]]
    print(f"Check failures: query length {bad_len}, X count {bad_x}, "
          f"templated with clip mismatches {bad_tmpl}, "
          f"splice_annotated below table identity {bad_ident}, "
          f"keys with read count != n_reads {len(off)}")
    if total != expected or off or bad_len or bad_x or bad_tmpl or bad_ident:
        sys.exit("Checks failed")

    for tmp, final in ((tmp_real, real_bam), (tmp_seg, seg_bam)):
        pysam.sort("-o", final, tmp)
        pysam.index(final)
        os.remove(tmp)
        print(f"Written: {final}")
    print(f"Written: {tsv}")


if __name__ == "__main__":
    main()
