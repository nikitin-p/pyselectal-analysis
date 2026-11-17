#!/usr/bin/env python
"""Stage 22 step 8: every long-clip read as STAR aligned it, with its soft clips drawn.

UCSC does not draw soft-clipped bases of a BAM record, so a long clip is invisible
in a BAM track of the original alignment. This step rewrites each record of
longclip_realigned.bam from its original STAR alignment (tags oc = CIGAR,
op = 1-based start) as one BED12 item:

  * every S run becomes a block on the bases directly beside the aligned part,
    where STAR left it, and is drawn thin (outside thickStart..thickEnd);
  * the aligned part (M, =, X, D) is drawn thick, and N gaps are drawn as lines;
  * the name is the original CIGAR, so the browser prints it beside the read in
    pack mode, and itemRgb is the read's origin-class colour (YC).

A clip that would run off the chromosome is cut at its end. The BED is sorted
and, when bedToBigBed is given, converted to bigBed.

Outputs in --outdir:
    longclip_cigar.bed   BED12, one item per read
    longclip_cigar.bb    the same as bigBed (only with --bedtobigbed)

Usage (from project root):
    python scripts/14_longclip/08_cigar_bed.py \\
        --realigned results/longclip/bam/longclip_realigned.bam \\
        --fai data/sacCer3.fa.fai --outdir results/longclip/bam \\
        [--bedtobigbed PATH] [--force]
"""
import argparse
import os
import re
import subprocess
import sys

import pysam

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--realigned", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--bedtobigbed", default=None)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def cigar_blocks(cigar, start):
    """Return (blocks, thick_start, thick_end) for a CIGAR placed at 0-based start.

    Blocks are [start, end) intervals in genome order. A leading S sits left of
    the aligned start and a trailing S right of the aligned end; contiguous runs
    are merged so that an S block and its neighbouring M block form one block.
    """
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise ValueError(f"Unparsable CIGAR {cigar}")
    lead = ops[0][0] if ops and ops[0][1] == "S" else 0
    trail = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "S" else 0
    pos, blocks, cur = start, [], None
    for n, op in ops:
        if op in "M=XD":
            if cur is None:
                cur = [pos, pos]
            pos += n
            cur[1] = pos
        elif op == "N":
            if cur is not None:
                blocks.append(cur)
                cur = None
            pos += n
    if cur is not None:
        blocks.append(cur)
    thick_start, thick_end = start, pos
    if lead:
        blocks[0][0] -= lead
    if trail:
        blocks[-1][1] += trail
    return blocks, thick_start, thick_end


def main():
    a = parse_args()
    bed_f = os.path.join(a.outdir, "longclip_cigar.bed")
    bb_f = os.path.join(a.outdir, "longclip_cigar.bb")
    target = bb_f if a.bedtobigbed else bed_f
    if os.path.exists(target) and not a.force:
        print(f"{target} exists; skipping (use --force)", file=sys.stderr)
        return

    sizes = {}
    with open(a.fai) as f:
        for line in f:
            chrom, size = line.split("\t")[:2]
            sizes[chrom] = int(size)

    rows, n_clipped = [], 0
    with pysam.AlignmentFile(a.realigned) as bam:
        for r in bam:
            cigar, op = r.get_tag("oc"), int(r.get_tag("op")) - 1
            chrom = r.reference_name
            blocks, ts, te = cigar_blocks(cigar, op)
            blocks[0][0] = max(blocks[0][0], 0)
            blocks[-1][1] = min(blocks[-1][1], sizes[chrom])
            n_clipped += "S" in cigar
            cs, ce = blocks[0][0], blocks[-1][1]
            rows.append((chrom, cs, ce, cigar, 0, "-" if r.is_reverse else "+", ts, te,
                         r.get_tag("YC"), len(blocks),
                         ",".join(str(e - s) for s, e in blocks) + ",",
                         ",".join(str(s - cs) for s, _ in blocks) + ","))
    if n_clipped != len(rows):
        raise SystemExit(f"{len(rows) - n_clipped} records carry no S in their original CIGAR")

    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(bed_f, "w") as f:
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"Written: {bed_f} ({len(rows)} reads)", file=sys.stderr)

    if a.bedtobigbed:
        chrom_sizes = os.path.join(a.outdir, "chrom.sizes.tmp")
        with open(chrom_sizes, "w") as f:
            for c, s in sizes.items():
                f.write(f"{c}\t{s}\n")
        subprocess.run([a.bedtobigbed, "-type=bed12", bed_f, chrom_sizes, bb_f], check=True)
        os.remove(chrom_sizes)
        print(f"Written: {bb_f}", file=sys.stderr)


if __name__ == "__main__":
    main()
