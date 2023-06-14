#!/usr/bin/env python3

import pandas as pd
import sys

library = sys.argv[1]
scalefactor = float(sys.argv[2])
outname = sys.argv[3]

def modify_value(line, scalefactor):
    # Split the line into a list of values using tabs as the separator    
    values = line.split()
    if len(values) >= 3:
        values[2] = str(int(values[2])*float(scalefactor))
    # Join the modified values back into a single line using tabs
    modified_line = "\t".join(values)
    return modified_line+"\n"


def scale_lib(library,scalefactor,filename):
    with open(filename, 'w') as file:
        file1 = open(library, 'r')
        Lines = file1.readlines()

        # INIT
        header_finished = False
        
        for line in Lines: 
            
            if line.startswith("#"):
                file.write(line)
                header_finished = True
                continue    
            if not header_finished: 
                file.write(line)
                continue             
            if line.startswith("!") and header_finished == True:
                file.write(line)
                continue
            else:
                if not line.startswith("#"): 
                    file.write(modify_value(line,scalefactor))
scale_lib(library, scalefactor, outname)