#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd

def get_seq_lengths(fasta_file):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    summary = pd.DataFrame()
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        l = len(sequence)
        name = name.replace("/", "_")
        entry = pd.DataFrame([{"family": "${fam_name}",
                               "sequence": name,
                               "sequence length": l,
                               "benchmarking_dataset": "${dataset}"}])
        summary = pd.concat([summary,entry], ignore_index=True)
    return(summary)

summary_lengths = get_seq_lengths("${fasta}")
summary_lengths.to_csv("${dataset}_${fam_name}_lengths.csv", index = False)
