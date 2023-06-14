#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys 
input_file = sys.argv[1]
output_file = sys.argv[2]
output_file_2 = sys.argv[3]


sequences = {}

with open(input_file) as f, open(output_file, 'w') as out, open(output_file_2, 'w') as out_nogap:    
    for line in f.readlines(): 
        if((line[0].isalpha()) & (~line.startswith("CLUSTAL") )) :
            seq_id = line.split()[0]
            seq = line.split()[1]
            if seq_id not in sequences:
                sequences[seq_id] = seq
            else:  
                sequences[seq_id] += seq
with open(output_file, 'w') as out, open(output_file_2, 'w') as out_nogap:    
    for seq_id in sequences.keys(): 
        seq = sequences[seq_id]
        seq_nogap = seq.replace("-", "")
        seq_id_clean = seq_id.replace("/", "_")
        seq_record = SeqRecord( seq=Seq(seq),
                                id=seq_id_clean,
                                description='')
        seq_record_nogap = SeqRecord( seq=Seq(seq_nogap),
                    id=seq_id_clean,
                    description='')
        SeqIO.write(seq_record, out, 'fasta')
        SeqIO.write(seq_record_nogap, out_nogap, 'fasta')           