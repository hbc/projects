#/usr/bin/env python
"""Prepare input configuration files for priority 1 variants.

Groups into family based calling for high priority variants.
"""
import collections
import csv
import datetime
import glob
import os
import sys

import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = _add_base_dir(yaml.load(in_handle))
    idremap = read_remap_file(config["idmapping"])
    baminfo = get_bam_files(config["inputs"], idremap)
    famsamples = priority_by_family(config["priority"], config["params"]["priority"])
    g1, g2 = split_families(famsamples, config["params"]["max_samples"])
    for name, g in [("g1", g1), ("g2", g2)]:
        write_config(g, baminfo, name, config["params"]["priority"], config["out"])

# ## Write output configurations

def write_config(g, baminfo, name, priority, outbase):
    outdir = os.path.dirname(outbase)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    base, ext = os.path.splitext(outbase)
    outfile = "%s-%s%s" % (base, name, ext)
    out = {"upload": {"dir": "../final"},
           "fc_date": datetime.datetime.now().strftime("%Y-%m-%d"),
           "fc_name": "alz-p%s-%s" % (priority, name),
           "details": []}
    for family, samples in g:
        for sample in samples:
            cur = {"algorithm": {"aligner": "bwa",
                                 "variantcaller": "gatk-haplotype",
                                 "quality_format": "Standard",
                                 "coverage_interval": "genome",
                                 "coverage_depth": "high"},
                   "analysis": "variant2",
                   "genome_build": "GRCh37",
                   "metadata": {"batch": str(family)},
                   "description": str(sample),
                   "files": [baminfo[sample]["bam"]]}
            out["details"].append(cur)
    with open(outfile, "w") as out_handle:
        yaml.dump(out, out_handle, allow_unicode=False, default_flow_style=False)
    return outfile

def split_families(famsamples, max_samples):
    """Split families into groups that fit for processing.
    """
    famsamples = sorted(famsamples.iteritems(), key=lambda xs: len(xs[1]), reverse=True)
    for f, s in famsamples:
        print f, s
    group1 = []
    group2 = []
    cur_size = 0
    for f, s in famsamples:
        if cur_size + len(s) <= max_samples:
            cur_size += len(s)
            group1.append((f, s))
        else:
            group2.append((f, s))
    for name, xs in [("total", famsamples), ("g1", group1), ("g2", group2)]:
        print name, sum([len(x) for _, x in xs])
    return group1, group2

# ## Priority file

def priority_by_family(in_file, target_priority):
    """Retrieve samples of a given priority, grouped by family.
    """
    samples = collections.defaultdict(list)
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        header = reader.next() # header
        for parts in reader:
            fam_id, sample_id, priority = parts[1:4]
            status_flag = parts[16]
            if status_flag != "Exclude" and priority and int(priority) == target_priority:
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

