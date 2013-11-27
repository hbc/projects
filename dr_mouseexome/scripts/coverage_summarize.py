#!/usr/bin/env python
"""Explore coverage summary to identify problematic samples and missed genes.
"""
import collections
import math
import operator
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main(in_file):
    df = pd.read_csv(in_file)
    #df = df.head(100)
    df = df[list(df.columns[:2]) + [x for x in df.columns if x.endswith("pct_nocov")]]
    #print(df.describe())
    plot_gene_coverages(df, "%s-regions.pdf" % os.path.splitext(in_file)[0])
    #df = df.apply(row_to_coverage, axis=1)
    plot_sample_coverages(df, "%s.pdf" % os.path.splitext(in_file)[0])
    #df = df.apply(normalize_row, axis=1)
    #plot_coverages(df, "%s-normalized.pdf" % os.path.splitext(in_file)[0])

def _half_i(n):
    if n % 2 == 0:
        return n // 2 - 1
    else:
        return n // 2

def plt_defaults():
    plt.rcParams['axes.grid'] = False
    plt.rcParams["xtick.major.size"] = 0
    plt.rcParams["xtick.minor.size"] = 0
    plt.rcParams["ytick.major.size"] = 0
    plt.rcParams["ytick.minor.size"] = 0
    plt.rcParams["font.family"] = "Verdana"
    plt.rcParams["font.monospace"] = "Courier"

def plot_sample_coverages(df, out_file):
    """Plot coverage of samples stratified by genes, identifying problem samples.
    """
    df = df[list(df.columns[2:])]
    df.columns = [x.replace("_pct_nocov", "") for x in df.columns]
    plt_defaults()
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams["axes.linewidth"] = 0
    plt.figure()
    num_cols = int(math.ceil(math.sqrt(len(df.columns))))
    while len(df.columns) % num_cols < int(math.ceil(num_cols / 2.0)):
        num_cols += 1
    num_rows = int(math.ceil(len(df.columns) / float(num_cols)))
    axes = df.hist(color="k", alpha=0.5,
                   layout=(num_rows, num_cols),
                   range=(0, 100), bins=100,
                   cumulative=True, histtype="step", normed=True,
                   grid=False,
                   ylabelsize=7, sharey=True,
                   sharex=True, xlabelsize=7, xrot=45, figsize=(11, 8))
    for i, xs in enumerate(axes):
        for j, x in enumerate(xs):
            if x.get_title():
                nc_area = area_under_nocoverage(df[x.get_title()])
                x.text(40, 0.35, nc_area, fontsize=12)
                x.set_ylim([0.0, 1.0])
                if i == _half_i(num_rows) and j == 0:
                    x.set_ylabel("Cumulative fraction of total genes")
    xs[_half_i(num_cols)].set_xlabel("Percentage of a gene uncovered by reads")
    plt.savefig(out_file)

def plot_gene_coverages(df, out_file):
    df.index = [x.split(";")[-1] for x in df.gene.values]
    df = df[list(df.columns[2:])].T.sort_index(axis=1, ascending=False)
    print df.ix[:5, :10]
    plt_defaults()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    df.boxplot(fontsize=7, grid=True, vert=False)
    ax.xaxis.grid(False)
    plt.xlabel("Percentage of uncovered bases (%s samples)" % len(df))
    plt.tight_layout()
    plt.savefig(out_file)

def area_under_nocoverage(col):
    """Calculate the area under a no-coverage curve as numerical metric of coverage.

    Assesses the overall coverage as a single number
    """
    hist, _ = np.histogram(col.values, normed=True, range=(0, 100), bins=100)
    cdf = np.cumsum(hist)
    return int(round(sum(cdf)))

def row_to_coverage(row):
    """Convert from no coverage to coverage.
    """
    update_vals = [100 - x for x in row.values[2:]]
    return pd.Series(list(row.values[:2]) + update_vals, index=list(row.index))

def normalize_row(row):
    row_median = np.median(row.values[2:])
    update_vals = [x - row_median for x in row.values[2:]]
    return pd.Series(list(row.values[:2]) + update_vals, index=list(row.index))

if __name__ == "__main__":
    main(*sys.argv[1:])
