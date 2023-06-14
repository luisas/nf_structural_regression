#!/usr/bin/env python3

import pandas as pd
import sys

library = sys.argv[1]
newval = int(sys.argv[2])
outname = sys.argv[3]

def modify_value(line, newval):
    # Split the line into a list of values using tabs as the separator    
    values = line.split()
    if len(values) >= 4:
        values[2] = str(newval)
        values[3] = str(newval)
    # Join the modified values back into a single line using tabs
    modified_line = "\t".join(values)
    return (modified_line+"\n")


def mod_lib(library, newval, filename): 
    with open(filename, 'w') as file:
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
            else:
                if not line.startswith("#"): 
                    line_print = modify_value(line, newval)
            file.write(line_print)

mod_lib(library, newval, outname)