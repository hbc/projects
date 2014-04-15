#!/usr/bin/env python
"""Prepare bcbio-nextgen configuration files for calling structural variation pilot.
"""
import csv
import glob
import os
import subprocess
import sys

from bcbio import utils

priority_file = "/n/hsphS10/hsphfs1/chb/projects/tanzi_ad/data/gwas/AD-Master-v2.csv"
bam_dir = "/n/hsphS10/hsphfs2/tanzi_recalled"
name = "alz-sv-pilot"

def main(family_file, region_file):
    families = get_families(family_file)
    samples = add_bams(get_samples(families, priority_file), bam_dir)
    config_file = write_config_file(samples, name, region_file)

def write_config_file(samples, name, region_file):
    meta_file = os.path.join(utils.safe_makedir(os.path.join(name, "config")),
                             "%s.csv" % name)
    bams = []
    with open(meta_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["samplename", "description", "batch", "sex",
                         "aligner", "mark_duplicates", "variantcaller", "svcaller", "variant_regions"])
        for sample in sorted(samples, key=lambda x: x["family"]):
            bams.append(sample["bam"])
            writer.writerow([os.path.basename(sample["bam"]), sample["sample"], sample["family"],
                             sample["gender"], "false", "false", "false", "lumpy;cn.mops", region_file])
    subprocess.check_call(["bcbio_nextgen.py", "-w", "template", "freebayes-variant",
                           meta_file] + bams)

def add_bams(samples, bam_dir):
    out = []
    for sample in samples:
        sample["bam"] = get_bam(sample["sample"], bam_dir)
        out.append(sample)
    return out

def get_bam(sample, bam_dir):
    bam_files = glob.glob(os.path.join(bam_dir, "*", "final", sample, "%s-*am" % sample))
    assert len(bam_files) > 0, "Did not find BAM files for %s: %s" % (sample, bam_files)
    if len(bam_files) > 1:
        bam_files = [x for x in bam_files if x.endswith(".bam")]
    return bam_files[0]

def get_samples(families, fname):
    samples = []
    with open(fname) as in_handle:
        reader = csv.reader(in_handle)
        reader.next()  # header
        for parts in reader:
            family_id, sample_id, priority = parts[1:4]
            status_flag = parts[16]
            if status_flag != "Exclude" and family_id in families:
                samples.append({"sample": sample_id, "family": family_id, "gender": _get_gender(parts[12])})
    return samples

def _get_gender(gender):
    if gender.lower() in ["m", "male", "1"]:
        return "male"
    elif gender.lower() in ["f", "female", "2"]:
        return "female"
    else:
        return ""

def get_families(in_file):
    families = set([])
    with open(in_file) as in_handle:
        for line in in_handle:
            families.add(line.strip())
    return families

if __name__ == "__main__":
    main(*sys.argv[1:])
