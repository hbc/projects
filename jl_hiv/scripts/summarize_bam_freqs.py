#!/usr/bin/env python
"""Print summary of raw base frequencies at each position in an input BAM file.

This is helpful for debugging purposes and to identify genome regions with high
variability.

Usage:
   summarize_bam_freqs.py <bam file> <count YAML file>
"""
import sys
import collections
from contextlib import closing

from Bio.Seq import Seq
import yaml
import pysam

def main(in_file, count_yaml):
    count_dict = get_seq_counts(count_yaml)
    with closing(pysam.Samfile(in_file, "rb")) as in_bam:
        for col in in_bam.pileup(stepper="all", max_depth=1000000000):
            read_freqs, depth = get_pileup_readfreqs(col, count_dict)
            print col.pos + 1, read_freqs, "depth", depth, "read-ends", percentage_read_ends(col)

def get_seq_counts(in_yaml):
    out = {}
    with open(in_yaml) as in_handle:
        for x in yaml.load(in_handle):
            out[x["seq"]] = x["count"]
    return out

def percentage_read_ends(col):
    total = 0
    ends = 0
    for read in col.pileups:
        total += 1
        if read.is_head or read.is_tail:
            ends += 1
    return "%0.1f" % (float(ends) / float(total) * 100.0)

def get_pileup_readfreqs(col, count_dict):
    counts = collections.defaultdict(int)
    for read in col.pileups:
        seq = read.alignment.seq
        if read.alignment.is_reverse:
            seq = str(Seq(seq).reverse_complement())
        counts[read.alignment.seq[read.qpos]] += count_dict[seq]
    counts = dict(counts)
    if counts.has_key("N"):
        del counts["N"]
    total = float(sum(counts.values()))
    out = []
    for base in sorted(counts.keys()):
        out.append([base, "%0.1f" % (float(counts[base]) * 100.0 / total)])
    return out, total

if __name__ == "__main__":
    main(*sys.argv[1:])
