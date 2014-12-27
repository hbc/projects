#/usr/bin/env python
"""Identify initial input BAM files for transfer.
"""
import collections
import csv
import datetime
import glob
import os
import sys
import subprocess

import yaml

def main(sample_csv, config_file, env):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
        config = _add_env_kvs(config, env)
        config = _add_base_dir(config)
    idremap = read_remap_file(config["idmapping"])
    baminfo = get_bam_files(config["inputs"], idremap)
    with open(sample_csv) as in_handle:
        in_handle.readline()  # header
        for sample_id in (l.split(",")[0] for l in in_handle):
            print baminfo.get(sample_id).get("bam")

# ## Directory and name remapping

def dir_to_sample(dname, idremap):
    vcf_file = os.path.join(dname, "Variations", "SNPs.vcf")
    bam_file = os.path.join(dname, "Assembly", "genome", "bam", "%s.bam" % os.path.split(dname)[-1])
    if not os.path.exists(bam_file):
        print "BAM file missing", bam_file
        return None
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                illumina_id = line.split("\t")[-1].replace("_POLY", "").rstrip()
                if idremap.get(illumina_id) is None:
                    print "Did not find remap", illumina_id
                return {"id": idremap.get(illumina_id), "dir": dname,
                        "illuminaid": illumina_id,
                        "bam": bam_file}
    raise ValueError("Did not find sample information in %s" % vcf_file)

def get_bam_files(fpats, idremap):
    out = {}
    for fpat in fpats:
        for dname in glob.glob(fpat):
            if os.path.isdir(dname):
                x = dir_to_sample(dname, idremap)
                if x:
                    out[x["id"]] = x
    return out

def read_remap_file(in_file):
    out = {}
    with open(in_file) as in_handle:
        in_handle.next() # header
        for line in in_handle:
            patient_id, illumina_id = line.rstrip().split()
            out[illumina_id] = patient_id
    return out

# ## Utils

def _add_env_kvs(config, env):
    """Add key values specific to the running environment.
    """
    for k, val in config[env].iteritems():
        config[k] = val
    return config

def _add_base_dir(config):
    for k in ["idmapping", "priority", "coverage", "fam"]:
        config[k] = os.path.join(config["base_dir"], config[k])
    return config

if __name__ == "__main__":
    main(*sys.argv[1:])

