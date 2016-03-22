#!/bin/env python
"""Subset multiple callers in an input set into combinations of calls.

Returns all single calls listed in callers plus all pairs.
"""
import itertools
import os
import sys

from pysam import VariantFile

from bcbio.variation import vcfutils

def main(in_file, callers):
    callers = callers.split(",")

    for caller in callers:
        subset_by_callers(in_file, [caller])
    for c1, c2 in itertools.combinations(callers, 2):
        subset_by_callers(in_file, [c1, c2])

def subset_by_callers(in_file, callers):
    out_file = "%s-%s.vcf" % (in_file.replace(".vcf", "").replace(".gz", ""), "_".join(callers))
    if not os.path.exists(out_file) and not os.path.exists(out_file + ".gz"):
        want_callers = set(callers)
        reader = VariantFile(in_file)
        writer = VariantFile(out_file, "w", header=reader.header)
        count = 0
        for rec in reader:
            cur_callers = set(rec.info["set"].split("-"))
            if len(cur_callers & want_callers) > 0:
                count += 1
                writer.write(rec)
        print callers, count
    return vcfutils.bgzip_and_index(out_file, {})

if __name__ == "__main__":
    main(*sys.argv[1:])
