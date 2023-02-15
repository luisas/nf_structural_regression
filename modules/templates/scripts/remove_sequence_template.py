#!/usr/bin/env python3

from Bio import SeqIO
import sys
import re

file =sys.argv[1]
seq_id = sys.argv[2].split(":")[0]

ff = open(file, "r")
lines = ff.readlines()

for line in lines:
    if bool(re.search('>'+seq_id+"\s_P_.*", line)):
        continue
    else:
        print(line, end = "")
