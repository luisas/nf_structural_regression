import pandas as pd
import sys

library = sys.argv[1]
aggfunc = sys.argv[2]
outname = sys.argv[3]

def aggsum(x):
    print(x)
    x = x.astype(float)
    x = x.sum()
    return(x.astype(str))

def aggmax(x):
    x = x.astype(float)
    x = x.max()
    return(str(x))

def print_chunk(entries, file, aggfunc):
    for i in entries: 
        if aggfunc == "sum": 
            df = pd.DataFrame(entries[i]).groupby([0,1]).agg({2: aggsum, 3: 'first', 4: 'first'}).reset_index()
        elif aggfunc == "max": 
            df = pd.DataFrame(entries[i]).groupby([0,1]).agg({2: aggmax, 3: 'first', 4: 'first'}).reset_index()
        for x in df.values.tolist(): 
            file.write("\t".join(x)+"\n")

def sum_lib(library,filename, aggfunc):
    with open(filename, 'w') as file:
        file1 = open(library, 'r')
        Lines = file1.readlines()

        
        # INIT
        header_finished = False
        entries = {}
        
        for line in Lines: 
            
            # Change of pair
            if line.startswith("#"):
                if header_finished:
                    print_chunk(entries, file,aggfunc) 
                    entries = {}
                file.write(line)
                header_finished = True
                continue    

            # Prefix 
            if not header_finished: 
                file.write(line)
                continue             
            # Suffix
            if line.startswith("!") and header_finished == True:
                print_chunk(entries, file, aggfunc) 
                entries = {}
                file.write(line)
                continue
            else:
                # Within each pair
                if not line.startswith("#"): 
                    residue = int(line.split()[0])
                    if residue in entries:
                        entries[residue].append(line.split())
                    else: 
                        entries[residue] = [line.split()]

sum_lib(library, outname, aggfunc)