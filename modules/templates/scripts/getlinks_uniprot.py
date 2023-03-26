import pandas as pd
import sys
import csv


template = sys.argv[1]
output = sys.argv[2]
ending =sys.argv[3]

df = pd.read_csv(template, sep='\s', header = None, engine = "python")
df["newname"] = df[0].str.replace(">", "").str.replace("/", "_")
df["line"] = df.apply(lambda row : "cp \"AF-"+row[2]+"-F1-model_v4."+ending+"\" \""+row["newname"]+"_ref."+ending+"\"", axis=1)
df["line"].to_csv(output, sep=" ", header=None, index=False, quoting=csv.QUOTE_NONE, escapechar=",")
