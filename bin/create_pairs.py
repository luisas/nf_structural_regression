#!/usr/bin/env python3

import sys
from Bio import SeqIO
from itertools import combinations
import pandas as pd 
file = sys.argv[1]
outfile = sys.argv[2]

ids = [r.id for r in SeqIO.parse(file, "fasta")]
df_ids = pd.DataFrame(list(combinations(ids, 2)))
df_ids.to_csv(outfile, index = False, header = False)