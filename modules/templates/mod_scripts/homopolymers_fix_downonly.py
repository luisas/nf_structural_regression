#!/usr/bin/env python3

import pandas as pd
import sys
from Bio import SeqIO
import glob
import numpy as np 


library = sys.argv[1]
windowsize = int(sys.argv[2])
outname = sys.argv[3]

alphabet = "ACDEFGHIKLMNPQRSTVWYX"

freq_aa = my_dict = {
    'A': 0.048503611971104234,
    'C': 0.03762704443818995,
    'D': 0.12771054195202802,
    'E': 0.02126528442317916,
    'F': 0.03314882571848516,
    'G': 0.03684523251086719,
    'H': 0.03867154517309316,
    'I': 0.020277074147043186,
    'K': 0.02105888607436595,
    'L': 0.08265315695656253,
    'M': 0.010445007349032117,
    'N': 0.03206054351565187,
    'P': 0.08381023860900022,
    'Q': 0.06208837602026457,
    'R': 0.04089814554210839,
    'S': 0.07782468649341714,
    'T': 0.01983925946774244,
    'V': 0.16124714638646528,
    'W': 0.02242236607561685,
    'Y': 0.021603027175782594,
    'X': 0.0
}

matrixm = [
    [6, -3, 1, 2, 3, -2, -2, -7, -3, -3, -10, -5, -1, 1, -4, -7, -5, -6, 0, -2, 0],
    [-3, 6, -2, -8, -5, -4, -4, -12, -13, 1, -14, 0, 0, 1, -1, 0, -8, 1, -7, -9, 0],
    [1, -2, 4, -3, 0, 1, 1, -3, -5, -4, -5, -2, 1, -1, -1, -4, -2, -3, -2, -2, 0],
    [2, -8, -3, 9, -2, -7, -4, -12, -10, -7, -17, -8, -6, -3, -8, -10, -10, -13, -6, -3, 0],
    [3, -5, 0, -2, 7, -3, -3, -5, 1, -3, -9, -5, -2, 2, -5, -8, -3, -7, 4, -4, 0],
    [-2, -4, 1, -7, -3, 6, 3, 0, -7, -7, -1, -2, -2, -4, 3, -3, 4, -6, -4, -2, 0],
    [-2, -4, 1, -4, -3, 3, 6, -4, -7, -6, -6, 0, -1, -3, 1, -3, -1, -5, -5, 3, 0],
    [-7, -12, -3, -12, -5, 0, -4, 8, -5, -11, 7, -7, -6, -6, -3, -9, 6, -12, -5, -8, 0],
    [-3, -13, -5, -10, 1, -7, -7, -5, 9, -11, -8, -12, -6, -5, -9, -14, -5, -15, 5, -8, 0],
    [-3, 1, -4, -7, -3, -7, -6, -11, -11, 6, -16, -3, -2, 2, -4, -4, -9, 0, -8, -9, 0],
    [-10, -14, -5, -17, -9, -1, -6, 7, -8, -16, 10, -9, -9, -10, -5, -10, 3, -16, -6, -9, 0],
    [-5, 0, -2, -8, -5, -2, 0, -7, -12, -3, -9, 7, 0, -2, 2, 3, -4, 0, -8, -5, 0],
    [-1, 0, 1, -6, -2, -2, -1, -6, -6, -2, -9, 0, 4, 0, 0, -2, -4, 0, -4, -5, 0],
    [1, 1, -1, -3, 2, -4, -3, -6, -5, 2, -10, -2, 0, 5, -2, -4, -5, -1, -2, -5, 0],
    [-4, -1, -1, -8, -5, 3, 1, -3, -9, -4, -5, 2, 0, -2, 6, 2, 0, -1, -6, -3, 0],
    [-7, 0, -4, -10, -8, -3, -3, -9, -14, -4, -10, 3, -2, -4, 2, 6, -6, 0, -11, -9, 0],
    [-5, -8, -2, -10, -3, 4, -1, 6, -5, -9, 3, -4, -4, -5, 0, -6, 8, -9, -5, -5, 0],
    [-6, 1, -3, -13, -7, -6, -5, -12, -15, 0, -16, 0, 0, -1, -1, 0, -9, 3, -10, -11, 0],
    [0, -7, -2, -6, 4, -4, -5, -5, 5, -8, -6, -8, -4, -2, -6, -11, -5, -10, 8, -6, 0],
    [-2, -9, -2, -3, -4, -2, 3, -8, -8, -9, -9, -5, -5, -5, -3, -9, -5, -11, -6, 9, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

columns = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']

# Create an empty dictionary to store the matrix
matrix = {}

# Populate the matrix dictionary
for col_idx, col_name in enumerate(columns):
    matrix[col_name] = {}
    for row_idx, row_name in enumerate(columns):
        matrix[col_name][row_name] = matrixm[col_idx][row_idx]


def get_exp_score(aa_i):
    expected_score = 0
    for aa in alphabet: 
        expected_score += freq_aa[aa] * matrix[aa_i][aa]
    return expected_score

def get_delta(seq, i):
    # i starts from 1
    i = i - 1 
    left_i = max(0,i-windowsize)
    right_i = min(len(seq), i+windowsize)
    local_avg = 0 
    range_j = range(left_i, right_i)
    for j in range_j: 
        if i !=j: 
            local_avg += matrix[seq[i]][seq[j]]
    local_avg = local_avg/len(range_j)
    exp_score = get_exp_score(seq[i])
    delta =  - local_avg * (1/2*windowsize) + exp_score
    return(delta)

def modify_value(line, newval):
    # Split the line into a list of values using tabs as the separator    
    values = line.split()
    if len(values) >= 3:
        values[2] = str(newval)
    # Join the modified values back into a single line using tabs
    modified_line = "\t".join(values)
    return modified_line+"\n"


with open(outname, 'w') as file:
    file1 = open(library, 'r')
    Lines = file1.readlines()
    count = 0
    header_finished = False

    seq_indexes = pd.DataFrame()
    residues = pd.DataFrame()   

    # Strips the newline character
    for line in Lines:
        line_print = line
        if line.startswith("!"):
            count += 1
            file.write(line_print)
            continue
        if count == 0:
            count += 1
            file.write(line_print)
            continue
        elif count == 1: 
            nseq = int(line)
            count += 1
            file.write(line_print)
            continue
        elif count == nseq+2: 
            header_finished = True

        if not header_finished: 
            count += 1
            # ----- PARSE INDEXES SEQUENCES 
            seq_id = line.split()[0]
            seq = line.split()[2]
            idx_seq = count-2
            entry = pd.DataFrame({"index":[idx_seq], "seq_name":[seq_id], "seq":[seq], "seq_3di": [seq]})
            seq_indexes = pd.concat([seq_indexes, entry])
        else:
            # ------- PARSE RESIDUES 
            if line.startswith("#"): 
                seq1 = line.replace("#","").split(" ")[0]
                seq2 = line.replace("#","").split(" ")[1]
            else:
                res1 = line.split()[0]
                res2 = line.split()[1]
                value = line.split()[2]
                #print("res1 = " + str(res1))
                #print("res2 = " + str(res2))
                seq_3di_1 = seq_indexes[seq_indexes["index"] == int(seq1)].seq_3di[0]
                seq_3di_2 = seq_indexes[seq_indexes["index"] == int(seq2)].seq_3di[0]
                res1_CHR = seq_3di_1[int(res1)-1]
                res2_CHR = seq_3di_2[int(res2)-1]
                #print(res1_CHR)
                #print(res2_CHR)
                #print("delta 1: " + str(get_delta(seq_3di_1, int(res1))))
                #print("delta 2: " + str(get_delta(seq_3di_2, int(res2))))
                s1 = max(0, int(res1)-windowsize)
                s2 = max(0, int(res2)-windowsize)
                e1 = min(len(seq_3di_1)-1, int(res1)+windowsize)
                e2 = min(len(seq_3di_2)-1, int(res2)+windowsize)
                #print(seq_3di_1[s1:e1])
                #print(seq_3di_2[s2:e2])
                delta = sum([get_delta(seq_3di_1, int(res1)), get_delta(seq_3di_2, int(res2))])/2
                #print("delta: " + str(delta))
                #print("substitution score "+ str(matrix[res1_CHR][res2_CHR]))
                perc_change = delta*100/(matrix[res1_CHR][res2_CHR]+0.1)* np.sign(matrix[res1_CHR][res2_CHR])
                #print(perc_change)
                if perc_change < 0: 
                    newval = max(0,int(value)*(1+perc_change/100))
                    newval = min(newval, 1000)
                else: 
                    newval = value 
                #print("current value: "+ str(value))
                #print("new val : "+ str(newval))
                #print("----")
                res_entry = pd.DataFrame({"seq1":[seq1], "seq2":[seq2], "res1":[res1], "res2":[res2], "value":[value]})
                residues = pd.concat([residues,res_entry])
                line_print = modify_value(line, newval)
        file.write(line_print)
