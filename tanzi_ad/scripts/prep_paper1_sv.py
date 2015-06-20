#!/usr/bin/env python
"""Prepare configuration file for running SV analysis on regions of interest for paper 1.

Usage:
  cd /n/hsphS10/hsphfs2/tanzi_recalled
  prep_paper1_sv.py <config_file>
"""
import csv
import math
import os
import subprocess
import sys

import toolz as tz
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(os.path.join(tz.get_in(["dirs", "inputs"], config),
                           tz.get_in(["inputs", "families"], config))) as in_handle:
        families = set([x.strip() for x in in_handle])
    fnames = read_base_config(tz.get_in(["resources", "basecfg"], config))
    samples = prep_sample_info(tz.get_in(["resources", "fam"], config), fnames, families)
    write_config_file(samples, tz.get_in(["dirs", "inputs"], config), tz.get_in(["inputs", "regions"], config),
                      tz.get_in(["dirs", "config"], config), tz.get_in(["configs", "template"], config))

def write_config_file(samples, input_dir, region_file, config_dir, template):
    region_file = os.path.relpath(_convert_to_bed(os.path.join(input_dir, region_file)), config_dir)
    template_file = os.path.join(config_dir, template)
    meta_file = os.path.join(config_dir, "%s.csv" % input_dir.split("/")[0])
    bams = []
    with open(meta_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        keys = ["description", "batch", "sex", "age", "apoe", "onset", "phenotype"]
        writer.writerow(["samplename"] + keys + ["variant_regions"])
        for sample in sorted(samples, key=lambda x: (x["batch"], x["description"])):
            bams.append(sample["file"])
            writer.writerow([os.path.basename(sample["file"])] + [sample[k] for k in keys] +
                            [region_file])
    subprocess.check_call(["bcbio_nextgen.py", "-w", "template", template_file, meta_file] + bams)

def prep_sample_info(fam_file, fnames, families):
    out = []
    with open(fam_file) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            fid, sid, _, _, sex, affected, _, _, age, onset, _, apoe, _ = line.split("\t")
            if fid in families and sid in fnames:
                out.append({"description": sid, "batch": fid, "file": fnames[sid],
                            "sex": _get_sex(sex),
                            "phenotype": _get_affected(affected),
                            "age": age, "onset": onset, "apoe": _get_apoe(apoe)})
    return out

def _convert_to_bed(in_file):
    """Convert file of regions to process into BED, rounding to 10kb on either side.
    """
    scale = 10000
    out_file = "%s.bed" % os.path.splitext(in_file)[0]
    with open(out_file, "w") as out_handle:
        with open(in_file) as in_handle:
            for line in in_handle:
                chrom, start, end = line.strip().split()
                start = int(math.floor(float(start) / scale) * scale)
                end = int(math.ceil(float(end) / scale) * scale)
                out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))
    return out_file

def _get_sex(sex):
    coding = {"1": "male", "2": "female"}
    return coding.get(str(sex), "unknown")

def _get_affected(affected):
    """0, -9: unknown; 1: unaffected 2: affected
    """
    coding = {"1": "unaffected", "2": "affected"}
    return coding.get(str(affected), "unknown")

def _get_apoe(x):
    coding = {"0": "not_e4", "1": "het_e4", "2": "hom_e4"}
    return coding.get(str(x), "unknown")

def read_base_config(in_file):
    with open(in_file) as in_handle:
        conf = yaml.load(in_handle)
    out = {}
    for sample in conf["details"]:
        out[sample["description"]] = sample["files"]
    return out

if __name__ == "__main__":
    main(sys.argv[1])
