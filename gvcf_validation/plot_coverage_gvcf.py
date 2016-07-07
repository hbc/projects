#!/usr/bin/env python
"""Prepare plots comparing coverage to gVCF based validations.
"""
import collections
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def inputs():
    Val = collections.namedtuple("Val", "callable,genome,fname")
    out = []
    for callable in ["coverage", "gvcf"]:
        for genome in ["NA12878", "NA24385"]:
            out.append(Val(callable, genome, "%s/grading-summary-%s.csv" % (callable, genome)))
    return out

def main():
    out_file = "gvcf_v_callable.pdf"
    df = plot_by_aligner_caller(out_file)
    out_file = "gvcf_v_callable-facet.pdf"
    plot_facet(df, out_file)

def plot_facet(df, out_file):
    with PdfPages(out_file) as pdf_out:
        for genome, genome_df in df.groupby(["genome"]):
            print genome_df.head()
            g = sns.FacetGrid(genome_df, col="aligner", row="caller", hue="best")
            g.map(sns.barplot, "difference", "metric")
            g.add_legend()
            pdf_out.savefig(g.fig)
            plt.clf()

def plot_by_aligner_caller(out_file):
    idfs = []
    for ival in inputs():
        df = pd.read_csv(ival.fname).rename(columns={"sample": "aligner"})
        df["callable"] = ival.callable
        df["genome"] = ival.genome
        idfs.append(df)
    df = pd.concat(idfs)
    sns.despine()
    sns.set(style="white")
    compare_dfs = []
    with PdfPages(out_file) as pdf_out:
        for (aligner, caller, genome), genome_df in df.groupby(["aligner", "caller", "genome"]):
            diff_data = {"metric": [], "difference": [], "best": []}
            for (vtype, metric), compare_df in genome_df.groupby(["vtype", "metric"]):
                coverage, gvcf = None, None
                for _, row in compare_df.iterrows():
                    if row["callable"] == "coverage":
                        coverage = row["value"]
                    elif row["callable"] == "gvcf":
                        gvcf = row["value"]
                assert gvcf and coverage, compare_df
                diff_data["metric"].append("%s %s" % (vtype, metric))
                diff_data["difference"].append(gvcf - coverage)
                if metric == "tp":
                    accurate = "gvcf" if gvcf > coverage else "coverage"
                else:
                    accurate = "gvcf" if gvcf < coverage else "coverage"
                diff_data["best"].append(accurate)
            plot_df = pd.DataFrame(diff_data)
            labels = sorted(list(plot_df["metric"].unique()))
            labels.reverse()
            plot_df["metric"].categories = labels
            g = sns.barplot(data=plot_df, x="difference", y="metric", hue="best")
            g.set_xlabel("difference (gVCF versus callable regions)")
            g.set_ylabel("")
            g.set_title("%s %s %s" % (genome, aligner, caller))
            g.figure.set_size_inches(6, 4)
            g.figure.tight_layout()
            pdf_out.savefig(g.figure)
            base, ext = os.path.splitext(out_file)
            cur_plot_out = "%s-%s-%s-%s.png" % (base, genome, aligner, caller)
            g.figure.savefig(cur_plot_out)
            plt.clf()
            plot_df["aligner"] = aligner
            plot_df["caller"] = caller
            plot_df["genome"] = genome
            compare_dfs.append(plot_df)
    return pd.concat(compare_dfs)


if __name__ == "__main__":
    main()
