#!/usr/bin/env python
"""Summarize statistics of call correctness across multiple parameters.

Usage:
  summarize_stats.py <stats file in YAML format>
"""
import sys

import yaml

from bcbio.summarize import print_summary_counts

def main(in_file):
    with open(in_file) as in_handle:
        stats = yaml.load(in_handle)
    regions = ["rt", "gag"]
    for region in regions:
        print "** %s" % region
        rstats = [(float(s["qual"]), float(s["kmer"]), s)
                  for s in stats if s["region"] == region]
        rstats.sort()
        for (_, _, info) in rstats:
            print_summary_counts(info)

if __name__ == "__main__":
    main(*sys.argv[1:])
