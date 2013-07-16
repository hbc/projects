#/usr/bin/env python
"""Prepare input configuration files for priority 1 variants.

Groups into family based calling for high priority variants.
"""
import collections
import csv
import glob
import os
import sys

import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = _add_base_dir(yaml.load(in_handle))
    idremap = read_remap_file(config["idmapping"])
    baminfo = get_bam_files(config["inputs"], idremap)
    famsamples = priority_by_family(config["priority"], 1)
    for family, samples in famsamples.iteritems():
        print family, samples
    print sum([len(x) for x in famsamples.itervalues()])

# ## Priority file

def priority_by_family(in_file, target_priority):
    """Retrieve samples of a given priority, grouped by family.
    """
    samples = collections.defaultdict(list)
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for _, fam_id, sample_id, priority in reader:
            if int(priority) == target_priority:
                samples[fam_id].append(sample_id)
    return dict(samples)

# ## Directory and name remapping

def dir_to_sample(dname, idremap):
    vcf_file = os.path.join(dname, "Variations", "SNPs.vcf")
    bam_file = os.path.join(dname, "Assembly", "genome", "bam", "%s.bam" % os.path.split(dname)[-1])
    assert os.path.exists(bam_file)
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                illumina_id = line.split("\t")[-1].replace("_POLY", "").rstrip()
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

def _add_base_dir(config):
    for k in ["idmapping", "priority", "out"]:
        config[k] = os.path.join(config["base_dir"], config[k])
    return config

if __name__ == "__main__":
    main(sys.argv[1])

