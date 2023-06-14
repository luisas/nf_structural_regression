#!/usr/bin/env python3

import pandas as pd
import sys

library1 = sys.argv[1]
library2 = sys.argv[2]
outname = sys.argv[3]

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
                entries[seq_pair]+= (line)
            else: 
                entries[seq_pair] = line
    return(entries)

def print_entries(entries_A, entries_B, file): 
    for i in entries_A: 
        file.write(entries_A[i])
        file.write(entries_B[i])  

def concat_libs(library_A, library_B, filename): 
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
                seq_pair = line
                print_entries(entries, entries_library_B, file )
                entries  ={}
                file.write(line)
                header_finished = True
                continue    
            # Prefix 
            if not header_finished:
                file.write(line)
                continue
            if line.startswith("!") and header_finished == True:
                print_entries(entries, entries_library_B, file )
                file.write(line)
                continue
            # Within each pair
            if not line.startswith("#"): 
                if seq_pair in entries:
                    entries[seq_pair]+= (line)
                else: 
                    entries[seq_pair] = line

concat_libs(library1, library2, outname) 
