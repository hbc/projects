"""Identify low complexity regions in GiaB reference materials.

Calculates percentage of calls at low GMS (< 50) regions in full calls
compared with array locations.

Requires gemini and pybedtools
"""
import os
import sys
import subprocess

from gemini import GeminiQuery
import pybedtools
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    for input in ["calls", "array", "orig_calls"]:
        identify_low_complexity(input, config[input], config["regions"])

def identify_low_complexity(name, in_vcf, in_bed):
    gms_thresh = 50.0
    subset_vcf = subset_by_region(name, in_vcf, in_bed)
    gemini_db = create_gemini_db(subset_vcf)
    print name
    gq = GeminiQuery(gemini_db)
    gq.run("SELECT count(*) from variants")
    total = list(gq)[0]["count(*)"]
    gq = GeminiQuery(gemini_db)
    gq.run("SELECT count(*) from variants WHERE gms_illumina < %s OR "
           "gms_solid < %s OR gms_iontorrent < %s" % (gms_thresh, gms_thresh, gms_thresh))
    low_gms = list(gq)[0]["count(*)"]
    print low_gms, total, "%.4f" % (float(low_gms) / float(total) * 100.0)

def create_gemini_db(in_vcf):
    out_file = "%s.db" % os.path.splitext(in_vcf)[0]
    if not os.path.exists(out_file):
        subprocess.check_call(["gemini", "load", "-v", in_vcf, out_file])
    return out_file

def subset_by_region(name, in_vcf, in_bed):
    out_file = "%s-subset.vcf" % name
    if not os.path.exists(out_file):
        tmp_file = out_file + ".callsonly"
        pybedtools.BedTool(in_vcf).intersect(in_bed).saveas(tmp_file)
        with open(out_file, "w") as out_handle:
            with open(in_vcf) as in_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        out_handle.write(line)
            with open(tmp_file) as in_handle:
                for line in in_handle:
                    out_handle.write(line)
        os.remove(tmp_file)
    return out_file

if __name__ == "__main__":
    main(sys.argv[1])
