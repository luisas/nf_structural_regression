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
import matplotlib as mpl


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

def cumsum(x): 
    return(list(x).cumsum())
def cum_var(x):
    ind_na = x.notna()
    nn = ind_na.cumsum()
    x[x.isna()] = 0
    cumsum(x^2) / (nn-1) - (cumsum(x))^2/(nn-1)/nn


def plot_scatter_perc(df1,df2,xlabel,ylabel,
                      palette = sns.dark_palette("#3399FF", reverse = True, as_cmap=True),
                      log = True, 
                      title = "regressive on homfam", hue_var = "n_sequences", metric = "tc", size_fig = 1): 
    sns.set_context("talk")
    f, ax = plt.subplots(figsize=(8*size_fig,6.4*size_fig ))
    
    # Prep df 
    df = df1.merge(df2, on = ["family","n_sequences"])
    if(log == True):
        df["n_sequences"] = np.log10( df["n_sequences"] )  
    
    hue = df["n_sequences"]
    
    # colorbar 
    norm = mpl.colors.Normalize( vmin=hue.min(), vmax=hue.max())
    sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
    sm.set_array([])

    metric_x = metric+"_x"
    metric_y = metric+"_y"
    # Plot 
    ax = sns.scatterplot(data = df, x = metric_x,
                    y = metric_y,
                    hue = hue_var,
                    s = 120,
                    palette = palette)
    

    # % above the line
    perc_y_better_than_x = (len(list(filter(lambda ele: ele == True, list(df[metric_x] <= df[metric_y])))) / len(list(df[metric_x]  >= df[metric_y] ))) * 100
    ax.get_legend().remove()
    
    # Color bar 
    cbar =ax.figure.colorbar(sm, ticks=np.log10([1000,10000,80_000]), format=mpl.ticker.ScalarFormatter())
    cbar.ax.set_yticklabels(["1K", "10K", "80K"]) 
    cbar.ax.set_ylabel('# of sequences  (log10 scale)', rotation=270, labelpad = 20, fontsize = "small")
    # Diagonal line
    ax.axline((1, 1), slope=1, ls="--", c=".2", lw = 0.8)
    
    # Axis limits
    plt.xlim([0, 100])
    plt.ylim([0, 100])

    
    # Axis labels
    ax.set(xlabel=xlabel,
           ylabel=ylabel,
           title = title + "\n metric: "+metric+"\n (n = "+str(len(df[metric_x] ))+") \n\n % y >= x  "+str(round(perc_y_better_than_x,1))+" \n")
    


