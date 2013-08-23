#!/usr/bin/env python
"""Explore coverage summary to identify problematic samples and missed genes.
"""
import os
import sys

import pandas as pd
import numpy
import matplotlib
import matplotlib.pyplot as plt

def main(in_file):
    df = pd.read_csv(in_file)
    #df = df.head(100)
    df = df[list(df.columns[:2]) + [x for x in df.columns if x.endswith("pct_nocov")]]
    print(df.describe())
    #df = df.apply(row_to_coverage, axis=1)
    plot_coverages(df, "%s.pdf" % os.path.splitext(in_file)[0])
    #df = df.apply(normalize_row, axis=1)
    #plot_coverages(df, "%s-normalized.pdf" % os.path.splitext(in_file)[0])

def plot_coverages(df, out_file):
    #matplotlib.rc_file(os.path.join(os.path.dirname(__file__), "ggplotish.rc"))
    df = df[list(df.columns[2:])]
    df.columns = [x.replace("_pct_nocov", "") for x in df.columns]
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['axes.grid'] = False
    plt.rcParams["xtick.major.size"] = 0
    plt.rcParams["xtick.minor.size"] = 0
    plt.rcParams["ytick.major.size"] = 0
    plt.rcParams["ytick.minor.size"] = 0
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.monospace"] = "Courier"
    plt.rcParams["axes.linewidth"] = 0
    #plt.rcParams["figure.facecolor"] = 0.75
    plt.figure()
    df.hist(color="k", alpha=0.5,
            bins=100, cumulative=True, histtype="step", normed=True,
            grid=False,
            ylabelsize=7, sharey=True,
            sharex=True, xlabelsize=7, xrot=45, figsize=(11, 8))
    plt.savefig(out_file)

def row_to_coverage(row):
    """Convert from no coverage to coverage.
    """
    update_vals = [100 - x for x in row.values[2:]]
    return pd.Series(list(row.values[:2]) + update_vals, index=list(row.index))

def normalize_row(row):
    row_median = numpy.median(row.values[2:])
    update_vals = [x - row_median for x in row.values[2:]]
    return pd.Series(list(row.values[:2]) + update_vals, index=list(row.index))

if __name__ == "__main__":
    main(*sys.argv[1:])
