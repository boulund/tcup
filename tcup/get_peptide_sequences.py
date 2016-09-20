#!/usr/bin/env python3.5

from sys import argv, exit
from utils import read_fasta as read_fasta


if len(argv) < 2:
    print("usage: script.py FASTA discriminative_peptides.txt")
    exit(1)

fastafile = argv[1]
discpeps = argv[2]


disc = {}
with open(discpeps) as f:
    for line in f:
        if line.startswith("-"):
            break
        peptide, families = line.strip().split(" ", maxsplit=1)
        translation = str.maketrans("{'}", "   ")
        families = families.translate(translation).split(",")
        disc[peptide] = [fam.strip() for fam in families]


disc_seqs = {}
for header, sequence in read_fasta(fastafile):
    try:
        disc_seqs[sequence] = ", ".join(disc[header])
    except KeyError:
        pass

for seq, target in disc_seqs.items():
    print(seq, target, sep="\t")
