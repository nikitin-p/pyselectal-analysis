#!/usr/bin/env python3
"""
Generate small synthetic BAMs with a known mix of 5' end types.
Used for pipeline testing when real data are not available.

Usage:
    python scripts/00_setup/make_test_bam.py

Outputs (in test/bam/):
  cond_A_rep1.bam  cond_A_rep2.bam   -- condition A, two replicates (similar composition)
  cond_B_rep1.bam  cond_B_rep2.bam   -- condition B, enriched for 2Sg
  cond_C_rep1.bam                    -- condition C, mostly mapped 5' ends

All reads map to a 10 kb synthetic chromosome "chrTest".
"""

import os
import random
import pysam

OUTDIR = "test/bam"
CHROM  = "chrTest"
CHROM_LEN = 10_000
READ_LEN  = 36
SEED = 42

HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
    "PG": [{"ID": "make_test_bam", "PN": "make_test_bam"}],
}

BASE_POOL = "ACGT"

def rand_seq(n, rng):
    return "".join(rng.choice(list(BASE_POOL)) for _ in range(n))


def make_alignment(header, qname, pos, cigar_str, seq, rng, is_reverse=False, mapq=255):
    a = pysam.AlignedSegment(header)
    a.query_name = qname
    a.query_sequence = seq
    a.flag = 0x10 if is_reverse else 0
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = mapq
    a.cigar = pysam.AlignedSegment.parse_cigar(cigar_str) if isinstance(cigar_str, str) else cigar_str
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


def cigar_tuples_to_str(tuples):
    ops = {0: "M", 4: "S"}
    return "".join(f"{l}{ops[op]}" for op, l in tuples)


def write_sample(path, type_counts, rng):
    """
    type_counts: dict mapping end-type label to count.
      Labels: '1Sg', '2Sgg', '2Sgc', 'M' (fully mapped 5' end)
    """
    records = []
    qnum = 0

    for label, n in type_counts.items():
        for _ in range(n):
            pos = rng.randint(100, CHROM_LEN - READ_LEN - 10)
            qname = f"r{qnum:06d}"
            qnum += 1

            if label == "M":
                seq = rand_seq(READ_LEN, rng)
                cigar = [(0, READ_LEN)]  # 36M
            elif label.startswith("1Sg"):
                clip_base = label[3] if len(label) > 3 else rng.choice(list("acgt"))
                seq = clip_base + rand_seq(READ_LEN - 1, rng)
                cigar = [(4, 1), (0, READ_LEN - 1)]  # 1S35M
            elif label.startswith("2Sg"):
                clip_seq = (label[3:5] if len(label) > 4 else
                            rng.choice(list("acgt")) + rng.choice(list("acgt")))
                seq = clip_seq + rand_seq(READ_LEN - 2, rng)
                cigar = [(4, 2), (0, READ_LEN - 2)]  # 2S34M
            elif label.startswith("3S"):
                clip_seq = rand_seq(3, rng)
                seq = clip_seq + rand_seq(READ_LEN - 3, rng)
                cigar = [(4, 3), (0, READ_LEN - 3)]
            else:
                seq = rand_seq(READ_LEN, rng)
                cigar = [(0, READ_LEN)]

            records.append((pos, qname, cigar, seq))

    records.sort(key=lambda x: x[0])

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        hdr = bam.header
        for pos, qname, cigar, seq in records:
            a = pysam.AlignedSegment(hdr)
            a.query_name = qname
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 255
            a.cigar = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam.write(a)

    pysam.sort("-o", path + ".tmp.bam", path)
    os.replace(path + ".tmp.bam", path)
    pysam.index(path)
    print(f"  wrote {sum(type_counts.values()):>5} reads → {path}")


def main():
    rng = random.Random(SEED)
    os.makedirs(OUTDIR, exist_ok=True)

    samples = {
        # Condition A: ~60 % 1Sg, ~20 % 2Sg, ~20 % M — two similar replicates
        "cond_A_rep1": {"1Sg": 600, "2Sgg": 110, "2Sgc": 90, "M": 200},
        "cond_A_rep2": {"1Sg": 580, "2Sgg": 120, "2Sgc": 80, "M": 220},
        # Condition B: shifted towards 2Sg
        "cond_B_rep1": {"1Sg": 350, "2Sgg": 300, "2Sgc": 150, "M": 200},
        "cond_B_rep2": {"1Sg": 340, "2Sgg": 310, "2Sgc": 160, "M": 190},
        # Condition C: mostly mapped ends
        "cond_C_rep1": {"1Sg": 150, "2Sgg":  50, "2Sgc":  50, "M": 750},
    }

    print("Generating synthetic BAMs …")
    for name, counts in samples.items():
        write_sample(f"{OUTDIR}/{name}.bam", counts, rng)
    print("Done.")


if __name__ == "__main__":
    main()
