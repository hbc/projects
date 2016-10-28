#!/usr/bin/env python
"""Summarize UMI tags by genome positions.

Uses UMI RX tags + position mapping to describe duplication rates.
Assumes a sorted input file.

Usage:
  summarize_umis.py <bam_file>
"""
from __future__ import print_function
import collections
import os
import sys

import numpy as np
import pysam

def main(bam_file):
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    mapped = 0
    unmapped = 0
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam_iter:
        for rec in bam_iter:
            umi = rec.get_tag("RX")
            if rec.is_unmapped:
                unmapped += 1
                key = (None, None)
            else:
                mapped += 1
                chrom = bam_iter.getrname(rec.reference_id)
                pos = rec.reference_start
                key = (chrom, pos)
            counts[key][umi] += 1
    print("Mapped: %s (%.1f%%)" % (mapped, float(mapped) / (mapped + unmapped) * 100.0))
    print("  Positions: %s" % (len(counts) - 1))
    print("  Barcodes per position")
    umi_reduction = []
    for key in sorted(counts.keys()):
        if key[0] is not None:
            total_seqs = sum(counts[key].values())
            umi_count = len(counts[key])
            if umi_count > 2 or total_seqs > 10:
                umi_reduction.append(float(total_seqs) / umi_count)
                print("   % 2s % 10s -- UMIs: % 3s, Total sequences: % 4s, max seqs/UMI: % 3s" % 
                      (key[0], key[1], umi_count, total_seqs, max(counts[key].values())))
    print ("  Reduction -- median: %.1fx; max: %.1fx" % (np.median(umi_reduction), max(umi_reduction)))

    print("Unmapped: %s (%.1f%%)" % (unmapped, float(unmapped) / (mapped + unmapped) * 100.0))
    print("  Barcodes: %s" % len(counts[(None, None)].keys()))

if __name__ == "__main__":
    main(*sys.argv[1:])
