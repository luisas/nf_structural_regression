#!/usr/bin/env python3

from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
header_set = sys.argv[2].split(":")[0]

for seq_record in ffile:
    if header_set == seq_record.name:
        continue
    else:
        print(seq_record.format("fasta"))
