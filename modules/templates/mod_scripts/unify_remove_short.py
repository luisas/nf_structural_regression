#!/usr/bin/env python3

import pandas as pd
import sys

library = sys.argv[1]
newval = int(sys.argv[2])
outname = sys.argv[3]

def modify_value(line, newval):
    # Split the line into a list of values using tabs as the separator    
    values = line.split()
    if len(values) >= 3:
        values[2] = str(newval)
    # Join the modified values back into a single line using tabs
    modified_line = "\t".join(values)
    return (modified_line+"\n")


def mod_lib(library, newval, N, filename):
    with open(filename, 'w') as file:
        file1 = open(library, 'r')
        Lines = file1.readlines()
        count = 0
        header_finished = False
        seq_indexes = pd.DataFrame()
        residues = pd.DataFrame()
        # Strips the newline character
        current_stretch = 1
        prev_residue = -1
        line_stretch = ""
        prev_line= "prev"

        for line in Lines:
            line_print = line
            if line.startswith("!"):
                count += 1
                if current_stretch > N:                    
                    file.write(line_stretch)
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
                    residue = int(line.split()[0])
                    if residue == prev_residue + 1:
                        current_stretch += 1
                        line_print = modify_value(line, newval)
                        line_stretch = line_stretch+line_print
                    else:
                        if current_stretch > N:                    
                            file.write(line_stretch)
                        current_stretch = 1
                        line_print = modify_value(line, newval)
                        line_stretch = line_print
                    prev_residue = residue
                    prev_line = line
            file.write(line_print)

                

mod_lib(library, newval, 3, outname)