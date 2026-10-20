#!/usr/bin/env python3
"""Write a per-base BED4 of the genome sequence, for a letter track in the hub.

Each genome position becomes one 1-bp item whose name is the base. With
--complement the name is the complementary base instead, so the row of letters
is the reverse strand written in genome order;
read right to left, it is the reverse complement, i.e. the - strand 5'->3'.

In the hub the track is drawn dense, white, with "labelOnFeature on", which
reproduces the look of UCSC's own Base Position letters. Both halves of that
matter and were checked on the UCSC renderer: pack mode refuses to put two
touching items on one row, so 1-bp items zig-zag over two rows, and a coloured
item is drawn as a filled box with the letter in white on top.

Chromosomes are written in byte order of their names and positions in
ascending order, which is the order bedToBigBed requires, so the output needs
no separate sort.

Usage:
    python3 scripts/99_additional/make_sequence_bed.py data/sacCer3.fa out.bed [--complement]
"""
import sys

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.upper())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        sys.exit(__doc__)
    fasta, out = args
    complement = "--complement" in sys.argv

    seqs = read_fasta(fasta)
    with open(out, "w") as fh:
        for chrom in sorted(seqs, key=lambda c: c.encode()):
            seq = seqs[chrom]
            for i, base in enumerate(seq):
                letter = COMPLEMENT.get(base, "N") if complement else base
                fh.write(f"{chrom}\t{i}\t{i + 1}\t{letter}\n")


if __name__ == "__main__":
    main()
