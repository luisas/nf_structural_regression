#!/usr/bin/env python3
import pandas as pd
from os import listdir
from os.path import isfile, join
import numpy as np

def get_seq_lengths_stats():
    files = [f for f in listdir(".") if f.endswith('csv') ]
    summary = pd.DataFrame()
    for file in files:
        print(file)
        df = pd.read_csv(file)
        stats_df = df.groupby(by=['family']).agg({'sequence length':['mean', 'median', "max"]})["sequence length"].reset_index()
        stats_df["benchmarking_dataset"] = df["benchmarking_dataset"][1]
        summary = pd.concat(list([summary,stats_df]), ignore_index = True)
    return(summary)

summary_lengths = get_seq_lengths_stats()
summary_lengths.to_csv("summary_lengths.csv", index = False)
