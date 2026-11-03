"""NAD branch step 4b: merge per-run CTSS tables into position x run matrices.

Writes two gzipped TSVs, one per count column (umis and reads):
  chr pos strand base_p1 <run1> <run2> ...
Columns follow manifest order; a position absent from a run counts 0.
Positions whose total over all runs is below --min-total are dropped.

Usage:
  python make_matrix.py --manifest priority1_runs.tsv --ctss-dir DIR
                        --outdir DIR [--min-total 2]
"""
import argparse
import gzip
import os
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--ctss-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-total", type=int, default=2)
    a = ap.parse_args()

    with open(a.manifest) as f:
        runs = [line.split("\t")[0] for line in f.read().splitlines()[1:] if line.strip()]

    for col in ("umis", "reads"):
        counts = defaultdict(lambda: [0] * len(runs))
        base = {}
        for j, run in enumerate(runs):
            path = os.path.join(a.ctss_dir, f"{run}.ctss.tsv.gz")
            with gzip.open(path, "rt") as f:
                hdr = f.readline().rstrip("\n").split("\t")
                ci = hdr.index(col)
                for line in f:
                    x = line.rstrip("\n").split("\t")
                    key = (x[0], int(x[1]), x[2])
                    counts[key][j] = int(x[ci])
                    base[key] = x[3]
        out = os.path.join(a.outdir, f"matrix_{col}.tsv.gz")
        kept = 0
        with gzip.open(out, "wt") as o:
            o.write("chr\tpos\tstrand\tbase_p1\t" + "\t".join(runs) + "\n")
            for key in sorted(counts):
                v = counts[key]
                if sum(v) < a.min_total:
                    continue
                kept += 1
                o.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{base[key]}\t" + "\t".join(map(str, v)) + "\n")
        print(f"{col}: {kept} of {len(counts)} positions kept -> {out}")


if __name__ == "__main__":
    main()
