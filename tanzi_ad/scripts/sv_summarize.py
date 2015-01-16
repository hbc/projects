#!/usr/bin/env python
"""Summarize structural variant calls across multiple samples.

Plot a heatmap of structural variant calls shared across multiple samples
"""
import collections
import csv
import glob
import json
import operator
import os
import sys

import toolz as tz
import yaml
import matplotlib.pyplot as plt
import pandas as pd
import pybedtools
import seaborn as sns

from bcbio import utils
from bcbio.structural import ensemble
from bcbio.variation import bedutils

# Threshold variables based on manual examination of coverage plots
#
# Need to be in this percentage of total samples to include in final plot
TOTAL_AFFECTED_PCT = 0.25
# Need to be in this percentage of affected samples within a family to move along
AFFECTED_PCT = 0.60
# Unaffected samples with risk below this threshold are treated as true controls
UNAFFECTED_THRESH = 0.3
# Minimum number of affected samples to include in final plot
MIN_AFFECTED = 5
# Maximum size of CNV to pass on. Avoid CNV calls in deep sequenced repeat regions
MAX_CNV_COUNT = 50
# Maximum standard structural variant size to include. Avoids long error prone calls
# better captured with CNVs
MAX_SV_SIZE = 10000  # bp
# Maximum size of calls to include
MAX_SIZE = 300  # kb
# Minimum size of calls to include
MIN_SIZE = 400  # bp

def main(config_file, *ensemble_files):
    sample_info = _read_samples(config_file)
    samples = _order_sample_info(sample_info)
    ensemble_files, to_exclude = _read_manual_exclude(ensemble_files)
    combo_bed, sample_files = _combine_calls(ensemble_files, sample_info)
    call_df_csv, positions, batch_counts = _create_matrix(combo_bed, samples, sample_info, sample_files,
                                                          to_exclude)
    _plot_heatmap(call_df_csv, samples, positions, sample_info, batch_counts)

def _plot_heatmap(call_csv, samples, positions, sample_info, batch_counts):
    def sample_sort(x):
        batch = sample_info[x]["batch"]
        return (-batch_counts.get(batch, 0), batch, x)
    out_file = "%s.png" % os.path.splitext(call_csv)[0]
    df = pd.read_csv(call_csv)
    sv_rect = df.pivot(index="position", columns="sample", values="caller_support")
    sv_rect = sv_rect.reindex_axis(positions, axis=0)
    sv_rect = sv_rect.reindex_axis(["%s: %s" % (sample_info[x]["batch"], x)
                                    for x in sorted(samples, key=sample_sort)],
                                   axis=1)
    fig = plt.figure(tight_layout=True)
    plt.title("Shared structural variant calls for affected and unaffected in regions of interest",
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
    fig.set_size_inches(20, 8)
    fig.savefig(out_file)

def _create_matrix(combo_bed, samples, sample_info, sample_files, to_exclude):
    out_file = "%s-df.csv" % os.path.splitext(combo_bed)[0]
    toplot_file = "%s-toplot.yaml" % os.path.splitext(combo_bed)[0]
    positions = []
    to_plot = []
    batch_counts = collections.defaultdict(int)
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["position", "sample", "caller_support"])
        items = []
        with open(combo_bed) as in_handle:
            for line in in_handle:
                chrom, start, end, sample_str = line.strip().split("\t")
                size = (int(end) - int(start)) / 1000.0
                position = "%s:%s-%s (%.1fkb)" % (chrom, start, end, size)
                cur_samples = collections.defaultdict(list)
                for cursample in sample_str.split(","):
                    val = (1 if sample_info[cursample]["phenotype"] == "affected"
                           else float(sample_info[cursample]["risk"]) - 1.0)
                    cur_samples[cursample] = val
                if not (chrom, int(start), int(end)) in to_exclude:
                    if _pass_sv(cur_samples, samples, size):
                        items.append((int(chrom), int(start), position, dict(cur_samples)))
                        plot_data = _get_plot_data(cur_samples, samples, sample_info, size)
                        if plot_data:
                            #_summarize_events(chrom, int(start), int(end), plot_data, sample_files)
                            plot_data.update({"chrom": chrom, "start": int(start), "end": int(end)})
                            to_plot.append(plot_data)
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
                    if cur_count != 0:
                        batch_counts[sample_info[sample]["batch"]] += 1
    with open(toplot_file, "w") as out_handle:
        yaml.safe_dump(to_plot, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file, positions, dict(batch_counts)

def _summarize_events(chrom, start, end, plot_data, sample_files):
    """Provide a summary of events found within a region and samples of interest.
    """
    print "* %s:%s-%s (%sbp)" % (chrom, start, end, end - start)
    for batch, samples in sorted(plot_data["batches"].items()):
        print " - %s" % batch
        for descr, status in sorted(samples.items()):
            print "   => %s %s" % (descr, status)
            locs = pybedtools.BedTool(sample_files[descr]).intersect(
                pybedtools.BedTool("%s\t%s\t%s\n" % (chrom, start, end), from_string=True), u=True)
            for f in locs:
                if f.stop - f.start < (end - start) * 10:
                    print "      %s:%s-%s (%sbp) %s" % (f.chrom, f.start, f.stop, f.stop - f.start, f.name)

def _get_plot_data(cur_samples, samples, sample_info, kb_size):
    """Retrieve structural variants and families to plot coverage to evaluate events.
    """
    cases = [x for x, y in cur_samples.items() if y > 0]
    controls = [x for x, y in cur_samples.items() if y < 0]
    #if (len(cases) + len(controls) > MIN_AFFECTED
    #      and float(len(cases)) / float(len(cases) + len(controls)) > AFFECTED_PCT
    #      and kb_size <= MAX_SIZE):
    if True:
        batches = set([])
        out = {"batches": {}}
        for sid in cases + controls:
            batch = sample_info[sid]["batch"]
            batches.add(batch)
        for batch in sorted(list(batches)):
            samples = {}
            for sample, val in sample_info.items():
                if val["batch"] == batch:
                    phenotype = val["phenotype"]
                    if phenotype.startswith("unaffected"):
                        phenotype = "%s (%.2f)" % (phenotype, float(val["risk"]))
                    samples[sample] = phenotype
            out["batches"][batch] = samples
        return out

def _pass_sv(cur_samples, samples, kb_size):
    """Determine if we should pass a structural variant for visualization.

    Passes:
      - Likely affected variants: found in at least three affected individuals and
        at least 60% of occurrences in affected.
      - Protective variants: primarily present in only unaffected
      - Common variants: found in more than 60% of total samples.
    """
    n_case = len([x for x in cur_samples.values() if x > 0])
    n_control = len([x for x in cur_samples.values() if x < 0])
    return ((((n_case + n_control) >= (len(samples) * TOTAL_AFFECTED_PCT)) or
             (n_case > MIN_AFFECTED) or
             (n_control > MIN_AFFECTED))
            and kb_size <= MAX_SIZE)

def _read_samples(in_file):
    with open(in_file) as in_handle:
        config = yaml.load(in_handle)
    out = {}
    for data in config["details"]:
        out[str(data["description"])] = {"phenotype": tz.get_in(["metadata", "phenotype"], data),
                                         "batch": tz.get_in(["metadata", "batch"], data),
                                         "risk": tz.get_in(["metadata", "risk"], data)}
    return out

def _order_sample_info(info):
    def sort_by_batch(x):
        return (info[x]["batch"], info[x]["phenotype"], x)
    return sorted(info.keys(), key=sort_by_batch)

def _combine_calls(in_files, sample_info):
    sample_files = {}
    sample_dir = utils.safe_makedir(os.path.join(os.getcwd(), "prepped"))

    files_by_batch = collections.defaultdict(list)
    for fname in in_files:
        descr = os.path.basename(fname).split("-")[0]
        sample_files[descr] = fname
        cur_file = os.path.join(sample_dir, os.path.basename(fname))
        files_by_batch[sample_info[descr]["batch"]].append(cur_file)
        if not utils.file_uptodate(cur_file, fname):
            with open(cur_file, "w") as out_handle:
                with open(fname) as in_handle:
                    for line in in_handle:
                        chrom, start, end, calls = line.strip().split("\t")[:4]
                        info = {"sample": descr, "phenotype": sample_info[descr]["phenotype"],
                                "risk": sample_info[descr]["risk"],
                                "calls": calls.split(","), "region": "%s:%s-%s" % (chrom, start, end)}
                        is_cnv = any([x.startswith("cnv") for x in info["calls"]])
                        if int(end) - int(start) >= MIN_SIZE and (is_cnv or int(end) - int(start) < MAX_SV_SIZE):
                            out_handle.write("\t".join([chrom, start, end, json.dumps(info)]) + "\n")
    samples_by_batch = collections.defaultdict(list)
    for s, x in sample_info.items():
        x["sample"] = s
        samples_by_batch[x["batch"]].append(x)
    ready_files = []
    for batch, fnames in files_by_batch.items():
        merged_file = _merge_by_batch(batch, fnames)
        filtered_file = _filter_batch(merged_file, samples_by_batch[batch])
        ready_files.append(filtered_file)
    return ensemble.combine_bed_by_size(ready_files, "combined", os.getcwd(), {}), sample_files

def _filter_batch(in_file, all_samples):
    """Filter a batch based on finding support in affected and not unaffected individuals.
    """
    out_file = "%s-filtered%s" % os.path.splitext(in_file)
    if not utils.file_uptodate(out_file, in_file):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for line in in_handle:
                    chrom, start, end, support = line.strip().split("\t")
                    samples = _check_support([json.loads(x) for x in support.split("&&")], all_samples)
                    if samples:
                        out_handle.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ",".join(samples)))
    return out_file

def _check_support(calls, all_samples):
    """Determine if a set of sample calls has enough evidence to pass along the calls.

    Treats CNVs differently from other calls, since these are relative differences
    and will also be potentially present in unaffected samples.
    """
    counts = {"cnv": {"affected": [], "unaffected": []},
              "sv": {"affected": [], "unaffected": []}}
    orig_s = collections.defaultdict(list)
    for sample in all_samples:
        orig_s[sample["phenotype"]].append(sample["sample"])
    orig_s = dict(orig_s)
    risks = {}
    call_info = {}
    for call in calls:
        cnv_count = len([x for x in call["calls"] if x.startswith("cnv")])
        if cnv_count > 0:
            counts["cnv"][call["phenotype"]].append(call["sample"])
        if len(call["calls"]) - cnv_count > 0:
            counts["sv"][call["phenotype"]].append(call["sample"])
        risks[call["sample"]] = float(call["risk"])
        call_info[call["sample"]] = call
    cnv_pass_samples = _check_support_cnv(counts["cnv"], orig_s, call_info, risks)
    sv_pass_samples = _check_support_sv(counts["sv"], orig_s, risks)
    return sorted(list(set(cnv_pass_samples).union(sv_pass_samples)))

def _check_support_cnv(cur, orig, call_info, risks):
    """For CNV support check for support in affected individuals.
    """
    pass_samples = collections.defaultdict(list)
    total = 0
    extra_samples = []
    for sample in cur["affected"] + cur["unaffected"]:
        for call in call_info[sample]["calls"]:
            if call.startswith("cnv"):
                count = int(call.split("_")[0].replace("cnv", ""))
                if count < MAX_CNV_COUNT:
                    if risks[sample] >= UNAFFECTED_THRESH:
                        total += 1
                        pass_samples[count].append(sample)
                    else:
                        extra_samples.append(sample)
    if ((total >= len(orig["affected"]) * AFFECTED_PCT and len(pass_samples) > 1)
          or len(extra_samples) > 1):
        return reduce(operator.add, pass_samples.values(), []) + extra_samples
    return []

def _check_support_sv(cur, orig, risks):
    """Check for support in non-CNV called structural variants.

    Ignore calls that are also present in unaffected individuals
    unlikely to be affected.
    """
    if len(cur["affected"]) >= len(orig["affected"]) * AFFECTED_PCT:
        in_unaffected = False
        for x in cur["unaffected"]:
            if risks[x] < UNAFFECTED_THRESH:
                in_unaffected = True
        if not in_unaffected:
            return cur["affected"] + cur["unaffected"]
        elif len(cur["affected"]) == 0:
            return cur["unaffected"]
    return []

def _merge_by_batch(batch, fnames):
    """Merge all calls in a family into a single callset.
    """
    merge_dir = utils.safe_makedir(os.path.join(os.getcwd(), "merged"))
    clean_dir = utils.safe_makedir(os.path.join(merge_dir, "clean"))
    merge_file = os.path.join(merge_dir, "%s-ensemble.bed" % batch)
    if not utils.file_uptodate(merge_file, fnames[0]):
        for fname in glob.glob(os.path.join(merge_dir, "%s-ensemble*" % batch)):
            os.remove(fname)
    ensemble.combine_bed_by_size(fnames, batch, merge_dir, {}, delim="&&")
    return bedutils.clean_file(merge_file, {}, bedprep_dir=clean_dir)

def _read_manual_exclude(fnames):
    """Potentially find a file with manual exclusions.
    """
    to_exclude = set([])
    if not fnames[0].endswith(".bed"):
        exclude_file = fnames[0]
        fnames = fnames[1:]
        with open(exclude_file) as in_handle:
            for line in in_handle:
                if line.strip():
                    chrom, coords = line.strip().split(":")
                    start, end = coords.split("-")
                    to_exclude.add((chrom, int(start), int(end)))
    return fnames, to_exclude

if __name__ == "__main__":
    main(*sys.argv[1:])
