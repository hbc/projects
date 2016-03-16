#!/bin/env python
"""Normalize a VCF file with multiple call approaches to a single set of passing calls.
"""
import collections
import gzip
import sys

def main(in_file):
    out_file = "%s-cleaned.vcf" % in_file.replace(".vcf.gz", "")
    caller_count = collections.defaultdict(int)
    with gzip.open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            tumor_is = None
            for line in in_handle:
                if line.startswith("##"):
                    if not line.startswith("##contig"):
                        out_handle.write(line)
                elif line.startswith("#"):
                    tumor_is = get_tumor_indexes(line)
                    parts = line.strip().split("\t")
                    out = parts[:9] + ["TUMOR"]
                    out_handle.write("\t".join(out) + "\n")
                else:
                    assert tumor_is
                    parts = line.strip().split("\t")
                    final_call, callers = extract_tumor_calls(parts, tumor_is)
                    if final_call:
                        for caller in callers:
                            caller_count[caller] += 1
                        #parts[0] = chrom_to_hg19(parts[0])
                        #parts[7] += ";CALLERS=%s" % ",".join(callers)
                        out = parts[:9] + [final_call]
                        out_handle.write("\t".join(out) + "\n")
    print dict(caller_count)

def chrom_to_hg19(chrom):
    if chrom == "MT":
        return "chrM"
    else:
        return "chr%s" % chrom

def extract_tumor_calls(parts, tumor_is):
    calls = []
    callers = []
    for i, caller in tumor_is:
        gt = parts[i].split(":")[0]
        if gt != "./.":
            calls.append(parts[i])
            callers.append(caller)
    if len(calls) > 0:
        return calls[0], callers
    else:
        return None, None

def get_tumor_indexes(line):
    tumors = ["TUMOR", "24385-12878-30-200"]
    parts = line.strip().split("\t")
    format_i = parts.index("FORMAT")
    out = []
    for i, p in enumerate(parts[format_i+1:]):
        base, caller = p.split(".")
        if base in tumors:
            out.append((format_i + i + 1, caller))
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])

