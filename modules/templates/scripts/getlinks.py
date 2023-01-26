import pandas as pd
import sys
import csv


template = sys.argv[1]
output = sys.argv[2]
ending =sys.argv[3]

df = pd.read_csv(template, sep='\s', header = None, engine = "python")
df["newname"] = df[0].str.replace(">", "")
df["line"] = df.apply(lambda row : "cp \""+row[2]+"."+ending+"\" \""+row["newname"]+"_ref."+ending+"\"", axis=1)
df["line"].to_csv(output, sep=" ", header=None, index=False, quoting=csv.QUOTE_NONE, escapechar=",")
