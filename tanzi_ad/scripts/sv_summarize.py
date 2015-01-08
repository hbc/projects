#!/usr/bin/env python
"""Summarize structural variant calls across multiple samples.

Plot a heatmap of structural variant calls shared across multiple samples
"""
import collections
import csv
import os
import sys

import toolz as tz
import yaml
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from bcbio import utils
from bcbio.structural import ensemble

def main(config_file, *ensemble_files):
    sample_info = _read_samples(config_file)
    samples = _order_sample_info(sample_info)
    combo_bed = _combine_calls(ensemble_files)
    call_df_csv, positions = _create_matrix(combo_bed, samples, sample_info)
    _plot_heatmap(call_df_csv, samples, positions, sample_info)

def _plot_heatmap(call_csv, samples, positions, sample_info):
    out_file = "%s.png" % os.path.splitext(call_csv)[0]
    df = pd.read_csv(call_csv)
    sv_rect = df.pivot(index="position", columns="sample", values="caller_support")
    sv_rect = sv_rect.reindex_axis(positions, axis=0)
    sv_rect = sv_rect.reindex_axis(["%s: %s" % (sample_info[x]["batch"], x) for x in samples], axis=1)
    fig = plt.figure(tight_layout=True)
    plt.title("Shared structural variant calls for affected and unaffected individuals in families",
              fontsize=16)
    ax = sns.heatmap(sv_rect, cbar=False,
                     cmap=sns.diverging_palette(255, 1, n=3, as_cmap=True))
    colors = sns.diverging_palette(255, 1, n=3)
    b1 = plt.bar(0, 0, bottom=-100, color=colors[-1])
    b2 = plt.bar(0, 0, bottom=-100, color=colors[0])
    ax.legend([b1, b2], ["affected", "unaffected"], ncol=2,
              bbox_to_anchor=(0.85, 0.995), loc=3)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    fig.set_size_inches(25, 18)
    fig.savefig(out_file)

def _create_matrix(combo_bed, samples, sample_info):
    out_file = "%s-df.csv" % os.path.splitext(combo_bed)[0]
    positions = []
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["position", "sample", "caller_support"])
        items = []
        with open(combo_bed) as in_handle:
            for line in in_handle:
                chrom, start, end, sample_str = line.strip().split("\t")
                position = "%s:%s-%s (%.1fkb)" % (chrom, start, end, (int(end) - int(start)) / 1000.0)
                cur_samples = collections.defaultdict(list)
                for cursample_str in sample_str.split(","):
                    cursample, count = cursample_str.split("_")
                    val = 1 if sample_info[cursample]["phenotype"] == "affected" else -1
                    cur_samples[cursample] = val
                if _pass_sv(cur_samples, samples):
                    items.append((int(chrom), int(start), position, dict(cur_samples)))
        items.sort()
        added = set([])
        for _, _, position, cur_samples in items:
            for i, sample in enumerate(samples):
                cur_count = cur_samples.get(str(sample), 0)
                key = (position, sample)
                if key not in added:
                    if i == 0:
                        positions.append(position)
                    added.add(key)
                    writer.writerow([position, "%s: %s" % (sample_info[sample]["batch"], sample), cur_count])
    return out_file, positions

def _pass_sv(cur_samples, samples):
    """Determine if we should pass a structural variant for visualization.

    Passes:
      - Likely affected variants: found in at least three affected individuals and
        at least 60% of occurrences in affected.
      - Common variants: found in more than 60% of total samples.
    """
    min_affected = 3
    pct_affected = 0.60
    n_case = len([x for x in cur_samples.values() if x > 0])
    n_control = len([x for x in cur_samples.values() if x < 0])
    return (((n_case + n_control) >= (len(samples) * pct_affected)) or
            (n_case > min_affected and float(n_case) / (n_case + n_control) > pct_affected))

def _read_samples(in_file):
    with open(in_file) as in_handle:
        config = yaml.load(in_handle)
    out = {}
    for data in config["details"]:
        out[str(data["description"])] = {"phenotype": tz.get_in(["metadata", "phenotype"], data),
                                         "batch": tz.get_in(["metadata", "batch"], data)}
    return out

def _order_sample_info(info):
    def sort_by_batch(x):
        return (info[x]["batch"], info[x]["phenotype"], x)
    return sorted(info.keys(), key=sort_by_batch)

def _combine_calls(in_files):
    min_size = 250
    sample_dir = utils.safe_makedir(os.path.join(os.getcwd(), "prepped"))
    ready_files = []
    for fname in in_files:
        descr = os.path.basename(fname).split("-")[0]
        cur_file = os.path.join(sample_dir, os.path.basename(fname))
        with open(cur_file, "w") as out_handle:
            with open(fname) as in_handle:
                for line in in_handle:
                    chrom, start, end, calls = line.strip().split("\t")[:4]
                    supports = set([x.split("_")[-1] for x in calls.split(",")])
                    if int(end) - int(start) >= min_size:
                        name = "%s_%s" % (descr, len(supports))
                        out_handle.write("\t".join([chrom, start, end, name]) + "\n")
        ready_files.append(cur_file)
    return ensemble._combine_bed_by_size(ready_files, "combined", os.getcwd(), {})

if __name__ == "__main__":
    main(*sys.argv[1:])
