#!/usr/bin/env python
"""Generate ROC curves for various pipeline options.

Usage:
  roc_plots.py <config file>
"""
import sys
import collections

import yaml
import rpy2.robjects as rpy

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)["roc_plot"]
    df_py = collections.defaultdict(list)
    for curve in config["curves"]:
        for fp, tp, small_tp in read_stat_data(curve["file"], config):
            df_py["true.positive"].append(tp)
            df_py["false.positive"].append(fp)
            df_py["true.positive.small"].append(small_tp)
            df_py["group"].append(curve["name"])
    plot_roc_curves(df_py, config["file"])

def plot_roc_curves(df_py, out_file):
    vector_types = {"group": rpy.StrVector}
    df = {}
    for k, vals in df_py.iteritems():
        vtype = vector_types.get(k, rpy.FloatVector)
        df[k] = vtype(vals)
    rpy.r.assign("exp.data", rpy.r['data.frame'](**df))
    rpy.r.assign("out.file", out_file)
    rpy.r('''
    library(ggplot2)
    p <- ggplot(exp.data, aes(x=false.positive, y=true.positive, group=group))
    p <- p + geom_line(aes(colour=group))
    ggsave(out.file, p, width=5, height=5)
    ''')

def read_stat_data(in_file, config):
    region = config["region"]
    qual = config["qual"]
    small_thresh = config["small_thresh"]
    with open(in_file) as in_handle:
        file_info = yaml.load(in_handle)
    for call in (c for c in file_info
                 if c["qual"] == qual and c["region"] == region):
        yield _sum_call_values(call["calls"], small_thresh)

def _sum_call_values(calls, small_thresh):
    tp_right = 0
    tp_total = 0
    tp_right_small = 0
    tp_total_small = 0
    for c in calls:
        if c["percent"] == 100.0:
            wrong = c.get("partial", 0) + c.get("wrong", 0)
            total = c["correct"] + wrong
            fp = float(wrong) / total
        else:
            right = c.get("correct", 0)
            total = right + c.get("partial", 0) + c.get("wrong", 0)
            if c["percent"] < float(small_thresh):
                tp_right_small += right
                tp_total_small += total
            else:
                tp_right += right
                tp_total += total
    return fp, float(tp_right) / tp_total, float(tp_right_small) / tp_total_small

if __name__ == "__main__":
    main(*sys.argv[1:])
