#!/usr/bin/env python
"""Convert UMI tagged BAM aligned files

Usage:

    bams_to_insertionsites.py <BAM file 1> <BAM file 2> ...
"""
from __future__ import print_function
import argparse
import csv
import os
import sys

import pysam
import Levenshtein

def main(bam_files, edit_distance, mapqual_thresh, no_umis):
    for bam_file in bam_files:
        prepare_count_output(bam_file, edit_distance, mapqual_thresh, no_umis)

def get_umi_groups(bam_file, edit_distance):
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
        if edit_distance == 0:
            grs.append([cur_umi])
        else:
            for g in grs:
                if any(Levenshtein.distance(cur_umi, w) <= edit_distance for w in g):
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

def prepare_count_output(bam_file, edit_distance, mapqual_thresh, no_umis):
    out_file = "%s-counts.txt" % os.path.splitext(bam_file)[0]
    umi_groups = {} if no_umis else get_umi_groups(bam_file, edit_distance)
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam_iter:
        seen = set([])
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for rec in bam_iter:
                if not rec.is_unmapped and rec.mapping_quality > mapqual_thresh:
                    umi = rec.get_tag("RX") if umi_groups else ""
                    if umi is not None:
                        norm_umi = umi_groups[umi] if umi_groups else ""
                        chrom = bam_iter.getrname(rec.reference_id)
                        pos = rec.reference_start
                        if rec.is_reverse:
                            pos += len(rec.query_alignment_sequence)
                        key = (chrom, pos, norm_umi)
                        if key not in seen:
                            seen.add(key)
                            writer.writerow([chrom, pos, norm_umi, rec.query_alignment_sequence])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BAMs to insertion sites, with UMI collapsing")
    parser.add_argument("--umi-edit-distance", default=1, type=int,
                        help="Allowed mismaches between UMIs to collapse")
    parser.add_argument("--mapqual-thresh", default=15, type=int,
                        help="Threshold for filtering multi-mapping reads")
    parser.add_argument("--no-umis", default=False, action="store_true",
                        help="Indicate this input lacks UMIs")
    parser.add_argument("bam_files", nargs="+")
    if len(sys.argv) == 1:
        parser.print_help()
    args = parser.parse_args()
    main(args.bam_files, args.umi_edit_distance, args.mapqual_thresh, args.no_umis)
