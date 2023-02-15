import pandas as pd
import sys
import csv
import re


file = sys.argv[1]
full_name = sys.argv[2]
template = sys.argv[3]

# Extract ids in template
#print("Reading template")
ids_in_template = []
ff = open(template, "r")
lines = ff.readlines()
for line in lines:
    result = re.search(r">(.*)\s_P_.*", line)
    ids_in_template.append(result.group(1))
#print(ids_in_template)

tcs_df_complete = pd.DataFrame()
file_opened = open(file, "r")
for line in file_opened:
    if ":" in line and "CPU" not in line:
        s_line = line.split(sep = ":")
        seq = s_line[0].strip()
        #print(seq)
        if seq != "cons" and seq in ids_in_template:
            entry = pd.DataFrame({'msa':[full_name], 'sequence': seq,"tcs": s_line[1].strip()})
            tcs_df_complete =  pd.concat([tcs_df_complete, entry], ignore_index = True)

#print(tcs_df_complete)
min_tcs_sequence = tcs_df_complete.sort_values("tcs")["sequence"].iloc[0]
min_tcs_value = float(tcs_df_complete.sort_values("tcs")["tcs"].iloc[0])

print(str(min_tcs_sequence)+":"+str(min_tcs_value))
