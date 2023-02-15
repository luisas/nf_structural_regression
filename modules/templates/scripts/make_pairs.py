#!/usr/bin/env python3
# coding: utf-8

# In[1]:


from Bio import SeqIO
import sys




family_name = sys.argv[1]
input_file = sys.argv[2]



fasta_sequences = SeqIO.parse(open(input_file),'fasta')

# Extract all possible pairs in
ids = []
seqs = {}
for fasta in fasta_sequences:
    ids.append(fasta.id)
    seqs[fasta.id] = str(fasta.seq)

pairs_ids = [(a, b) for idx, a in enumerate(ids) for b in ids[idx + 1:]]


# Create new file
for pair in (pairs_ids):
    newname = family_name+"_pair_"+pair[0]+"_VS_"+pair[1]
    f = open(newname+".fa", "w")
    f.write(">"+pair[0] + '\n')
    f.write(seqs[pair[0]]+ '\n')

    f.write(">"+pair[1] + '\n')
    f.write(seqs[pair[1]]+ '\n')

    f.close()
