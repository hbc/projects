#!/usr/bin/env python
"""Prepare Illumina called variants, merging into single sample GATK-compatible VCFs.

Usage:
  illumina_variant_prep.py <in config> <cores>
"""
import csv
import glob
import multiprocessing
import pprint
import os
import subprocess
import sys

import yaml

import joblib

from bcbio.variation import effects, population, vcfutils
from bcbio.distributed.messaging import parallel_runner

def main(config_file, cores):
    cores = int(cores)
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    idremap = read_remap_file(config["idmapping"])
    exclude = read_priority_file(config["array"]["priority"], idremap)
    samples = list(get_input_samples(config["inputs"], idremap))
    problem = [x for x in samples if x["id"] is None]
    if len(problem) > 0:
        print "Problem identifiers"
        for p in problem:
            print p["illuminaid"], os.path.basename(p["dir"])
        raise NotImplementedError
    check_fam(samples, config["array"]["fam"])

    config["algorithm"] = {"num_cores": cores}
    out_files = [outf for outf in joblib.Parallel(cores)(joblib.delayed(run_illumina_prep)(s, config)
                                                         for s in samples if s["id"] is not None
                                                         and s["id"] not in exclude)]
    merge_file = merge_vcf_files(out_files, cores, config)
    effects_file = effects.snpeff_effects({"vrn_file": merge_file,
                                           "genome_build": "GRCh37",
                                           "config": config})
    noexclude_file = "%s-noexclude%s" % os.path.splitext(effects_file)
    noexclude_file = vcfutils.exclude_samples(effects_file, noexclude_file, exclude,
                                              config["ref"]["GRCh37"], config)
    gemini_db = population.prep_gemini_db([noexclude_file], "casava",
                                          [{"config": config, "work_bam": "yes", "genome_build": "GRCh37"}],
                                          {"config": config})[0][1]["db"]
    print gemini_db

def merge_vcf_files(sample_files, cores, config):
    out_file = config["outputs"]["merge"]
    config["algorithm"] = {}
    run_parallel = parallel_runner({"type": "local", "cores": min(cores, 8)}, {}, config)
    vcfutils.parallel_combine_variants(sample_files, out_file, config["ref"]["GRCh37"],
                                       config, run_parallel)
    return out_file

def check_fam(samples, fam_file):
    """Ensure identifiers are present in PLINK fam file.
    """
    fam_samples = set([])
    with open(fam_file) as in_handle:
        for line in in_handle:
            fam_samples.add(line.split()[1].strip())
    missing_ids = []
    for sample in samples:
        if sample["id"] not in fam_samples:
            missing_ids.append(sample["id"])
    with open("missing_fam_sample_ids.txt", "w") as out_handle:
        for x in sorted(missing_ids):
            out_handle.write("%s\n" % x)

def read_priority_file(in_file, idremap):
    """Read priority and failed information from input files.
    """
    exclude = []
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        header = reader.next()
        #for i, h in enumerate(header):
        #    print i, h
        for parts in reader:
            rutgers_id = parts[0]
            subject_id = parts[2]
            status_flag = parts[16]
            remap_s_id = idremap.get(rutgers_id)
            if status_flag == "Exclude":
                # if remap_s_id != subject_id:
                #     print rutgers_id, remap_s_id, subject_id
                #     print parts
                exclude.append(subject_id)
            else:
                assert status_flag in ["Use", "QC"], status_flag
    return set(exclude)

def run_illumina_prep(sample, config):
    tmp_dir = config.get("tmpdir", os.getcwd())
    if not os.path.exists(tmp_dir):
        try:
            os.makedirs(tmp_dir)
        except OSError:
            assert os.path.exists(tmp_dir)
    out_file = os.path.join(os.getcwd(), "%s.vcf" % sample["id"])
    if not os.path.exists(out_file):
        print sample["id"], sample["dir"], out_file
        subprocess.check_call(["java", "-Xms1g", "-Xmx2g", "-jar", config["bcbio.variation"],
                               "variant-utils", "illumina", sample["dir"],
                               sample["id"], config["ref"]["GRCh37"],
                               config["ref"]["hg19"],
                               "--outdir", os.getcwd(),
                               "--tmpdir", tmp_dir])
    return out_file

def dir_to_sample(dname, idremap):
    vcf_file = os.path.join(dname, "Variations", "SNPs.vcf")
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                illumina_id = line.split("\t")[-1].replace("_POLY", "").rstrip()
                return {"id": idremap.get(illumina_id), "dir": dname,
                        "illuminaid": illumina_id}
    raise ValueError("Did not find sample information in %s" % vcf_file)

def get_input_samples(fpats, idremap):
    for fpat in fpats:
        for dname in glob.glob(fpat):
            if os.path.isdir(dname):
                yield dir_to_sample(dname, idremap)

def read_remap_file(in_file):
    out = {}
    with open(in_file) as in_handle:
        in_handle.next() # header
        for line in in_handle:
            patient_id, illumina_id = line.rstrip().split()
            out[illumina_id] = patient_id
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
