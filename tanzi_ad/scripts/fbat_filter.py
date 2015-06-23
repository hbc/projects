#!/usr/bin/env python
"""Provide a filter for raw FreeBayes results as input to FBAT association testing.

Simple depth based filtering:

    - Keep only SNPs
    - Remove calls with low support across the sample
    - Remove low and high depth calls

Usage:
    fbat_filter.py <input_vcf> <out_file>
"""
import sys

from pysam import VariantFile

def main(in_file, out_file):
    config = {"callrate": 0.75,
              "min_depth": 13,
              "max_depth": 250}
    with VariantFile(in_file) as bcf_in:
        with VariantFile(out_file, "w", header=bcf_in.header) as bcf_out:
            for rec in bcf_in:
                if _passes(rec, config):
                    bcf_out.write(rec)

def _passes(rec, config):
    total = 0
    passing = 0
    for name, sample in rec.samples.items():
        total += 1
        depth = sample.get("DP", 0)
        if depth >= config["min_depth"] and depth < config["max_depth"]:
            passing += 1
    return float(passing) / float(total) > config["callrate"]


if __name__ == "__main__":
    main(*sys.argv[1:])
