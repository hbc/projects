#!/usr/bin/env python
"""Convert UMI tagged BAM aligned files

Usage:

    bams_to_insertionsites.py <BAM file 1> <BAM file 2> ...
"""
from __future__ import print_function
import csv
import os
import sys

import pysam

def main(*bam_files):
    for bam_file in bam_files:
        prepare_count_output(bam_file)

def prepare_count_output(bam_file):
    out_file = "%s-counts.txt" % os.path.splitext(bam_file)[0]
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam_iter:
        seen = set([])
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for rec in bam_iter:
                if not rec.is_unmapped:
                    umi = rec.get_tag("RX")
                    chrom = bam_iter.getrname(rec.reference_id)
                    pos = rec.reference_start
                    if rec.is_reverse:
                        pos += len(rec.query_alignment_sequence)
                    key = (chrom, pos, umi)
                    if key not in seen:
                        seen.add(key)
                        writer.writerow([chrom, pos, rec.query_alignment_sequence])

if __name__ == "__main__":
    main(*sys.argv[1:])
