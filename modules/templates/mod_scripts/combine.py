#!/usr/bin/env python3

import pandas as pd
import sys

library1 = sys.argv[1]
library2 = sys.argv[2]
aggfunc = sys.argv[3]
outname = sys.argv[4]


def aggsum(x):
    x = x.astype(float)
    x = x.sum()
    return(x.astype(str))

def aggmax(x):
    x = x.astype(float)
    x = x.max()
    return(str(x))



def lib_to_dict(library):
    file1 = open(library, 'r')
    Lines = file1.readlines()


    # INIT
    header_finished = False
    entries = {}

    for line in Lines: 

        # Change of pair
        if line.startswith("#"):
            seq_pair = line
            header_finished = True
            continue    
        # Prefix 
        if not header_finished: 
            continue
        if line.startswith("!") and header_finished == True:
            continue
        # Within each pair
        if not line.startswith("#"): 
            if seq_pair in entries:
                entries[seq_pair].append(line.split())
            else: 
                entries[seq_pair] =  [line.split()]
    return(entries)



def print_chunk(entries_A, entries_B, file, aggfunc = "sum"):
    entries = pd.concat([pd.DataFrame(entries_A), pd.DataFrame(entries_B)])
    if aggfunc == "sum": 
        df = entries.groupby([0,1]).agg({2: aggsum, 3: 'first', 4: 'first'}).reset_index()
    elif aggfunc == "max": 
        df = entries.groupby([0,1]).agg({2: aggmax, 3: 'first', 4: 'first'}).reset_index()
    # Sort 
    df[0] = df[0].astype(int)
    df[1] = df[1].astype(int)
    df = df.sort_values([0,1])
    df[0] = df[0].astype(str)
    df[1] = df[1].astype(str)   
    for x in df.values.tolist(): 
        file.write("\t".join(x)+"\n")


def print_entries(entries_A, entries_B, file): 
    for i in entries_A: 
        file.write(entries_A[i])
        file.write(entries_B[i])  

def combine_libs(library_A, library_B, aggfunc, filename): 
    with open(filename, 'w') as file:

        entries_library_B = lib_to_dict(library_B)

        file1 = open(library_A, 'r')
        Lines = file1.readlines()

        # INIT
        header_finished = False
        entries = {}

        for line in Lines: 

            # Change of pair
            if line.startswith("#"):
                if header_finished: 
                    print_chunk(entries[seq_pair], entries_library_B[seq_pair], file, aggfunc)
                entries  ={}
                seq_pair = line
                file.write(line)
                header_finished = True
                continue    
            # Prefix 
            if not header_finished:
                file.write(line)
                continue
            if line.startswith("!") and header_finished == True:
                print_chunk(entries[seq_pair], entries_library_B[seq_pair], file, aggfunc)
                file.write(line)
                continue
            # Within each pair
            if not line.startswith("#"): 
                if seq_pair in entries:
                    entries[seq_pair].append(line.split())
                else: 
                    entries[seq_pair] = [line.split()]



combine_libs(library1, library2, aggfunc, outname)