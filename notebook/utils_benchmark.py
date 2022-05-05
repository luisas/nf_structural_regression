def test():
    print("Import working!")

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import matplotlib.ticker as ticker
import sys


### parse trace files
def parse_trace(result_dir):
    f = os.path.basename(result_dir)
    name_trace = ".command.trace"
    trace = join(result_dir, name_trace)
    df = pd.read_csv(trace, sep = "=").T
    df = df.rename(index={'nextflow.trace/v2': f})
    return(df)

def split_name(df):
    df[["family", "method", "bucket_size", "dynamicX_label","dynamicX_val",
               "master_msa", "master_batch",
               "slave_msa", "slave_batch","tree"]] = df.name.str.split(".", expand = True)
    return(df)


### Collect all computation files across directories
def get_computation_times(evaluation_dir, dataset, task, extralevel = False, splitname =False):
    # Extract trace files w/ corresponding alignments
    traces = pd.DataFrame()
    alignments_dir = join(evaluation_dir,task)
    for fam in listdir(alignments_dir):
        family_dir=join(alignments_dir, fam)

        # remove extra level if needed
        if extralevel == False:
            family_dir = alignments_dir

        for f in listdir(family_dir):
            result_dir = join(family_dir,f)
            if(os.path.isdir(result_dir)):
                df = parse_trace(result_dir)
                traces = traces.append(df)

        if extralevel == False:
            break
    # Trace files parsed
    traces = traces.reset_index(level=0)
    traces = traces.rename(columns={'index': 'name'})
    traces["benchmarking_dataset"] = dataset
    if splitname:
        traces = split_name(traces)

    traces["task"] = task
    return(traces)



def get_evaluation(my_dir, csv_file, splitname = True):

    scores = join(my_dir, csv_file)
    scores_df = pd.read_csv(scores, sep=";", header = None).iloc[: , :-1]
    scores_df.set_axis(list(["name", "sp","tc", "column"]), axis=1, inplace=True)

    method = csv_file.split(".")[0]


    if splitname:
        if method == "dynamic":
            scores_df[["family", "method", "bucket_size", "dynamicX_label","dynamicX_val",
                       "master_msa",
                       "slave_msa","tree"]] = scores_df.name.str.split(".", expand = True)
        elif method == "regressive":
            scores_df[["family", "method", "bucket_size",
                       "master_msa","tree"]] = scores_df.name.str.split(".", expand = True)
        elif method == "progressive":
            scores_df[["family", "method",
                       "master_msa","tree"]] = scores_df.name.str.split(".", expand = True)
    #scores_df.master_msa = scores_df.master_msa.str.split("_")[0]
    #scores_df.master_msa = scores_df.master_msa.apply(lambda x: x.astype(str).str.upper())
    scores_df = scores_df.reset_index()
    return(scores_df)

def get_evaluation_all(my_dir):
    all_files = os.listdir(my_dir)
    csv_files = list(filter(lambda f: f.endswith('.csv'), all_files))
    df  = pd.DataFrame()
    for file in csv_files:
        df = pd.concat(list([df, get_evaluation(my_dir, file)]))
    return(df)

def cumavg(x): 
    return(x.cumsum()/list(range(len(x)+1)[1:]))
