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
import Levenshtein

MAPQUAL_THRESH = 15
EDIT_DISTANCE = 1

def main(*bam_files):
    for bam_file in bam_files:
        prepare_count_output(bam_file)

def get_umi_groups(bam_file):
    """Retrieve UMI groups with UMIs within a certain edit distance.

    Based on logic from:
    http://stackoverflow.com/a/35173198/252589
    """
    all_umis = set([])
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam_iter:
        for rec in bam_iter:
            all_umis.add(rec.get_tag("RX"))
    print(len(all_umis))
    grs = []
    for i, cur_umi in enumerate(sorted(all_umis)):
        if i % 1000 == 0:
            print(i, len(grs))
        for g in grs:
            if any(Levenshtein.distance(cur_umi, w) <= EDIT_DISTANCE for w in g):
                g.append(cur_umi)
                break
        else:
            grs.append([cur_umi])
    out = {}
    for cur_gr in grs:
        base = Levenshtein.median(cur_gr)
        for gr in cur_gr:
            out[gr] = base
    return out

def prepare_count_output(bam_file):
    out_file = "%s-counts.txt" % os.path.splitext(bam_file)[0]
    umi_groups = get_umi_groups(bam_file)
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam_iter:
        seen = set([])
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for rec in bam_iter:
                if not rec.is_unmapped and rec.mapping_quality > MAPQUAL_THRESH:
                    umi = rec.get_tag("RX")
                    if umi:
                        norm_umi = umi_groups[umi]
                        chrom = bam_iter.getrname(rec.reference_id)
                        pos = rec.reference_start
                        if rec.is_reverse:
                            pos += len(rec.query_alignment_sequence)
                        key = (chrom, pos, norm_umi)
                        if key not in seen:
                            seen.add(key)
                            writer.writerow([chrom, pos, norm_umi, rec.query_alignment_sequence])

if __name__ == "__main__":
    main(*sys.argv[1:])
