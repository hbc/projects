#!/usr/bin/env python
"""Prepare clean aligned files (only HIV data, no human) for submission to SRA.
"""
import os
import subprocess
import sys

import joblib
import yaml

def main(config_file, cores):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    out_dir = os.path.join(os.getcwd(), "sra")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    joblib.Parallel(cores)(joblib.delayed(generate_aligned_only)
                           (exp, config["name"], out_dir) for exp in config["experiments"])

def generate_aligned_only(exp, exp_name, out_dir):
    if exp["description"] == "Control":
        sample_name = exp["description"]
    else:
        sample_name = os.path.basename(exp["files"][0]).split("_")[0]
    base_out = os.path.join(out_dir, "%s-%s" % (exp_name, sample_name))

    fq1 = "%s-1.fq" % base_out
    fq2 = "%s-2.fq" % base_out
    fastq_file, pair_file = exp["files"]
    ref_file = exp["ref"]
    cmd = ("novoalign -o SAM -d {ref_file} -f {fastq_file} {pair_file} "
           "| samtools view -b -S -u - "
           "| bamtofastq F={fq1} F2={fq2} S=/dev/null O=/dev/null O2=/dev/null")
    if not os.path.exists(fq1):
        subprocess.check_call(cmd.format(**locals()), shell=True)
    for fq in [fq1, fq2]:
        if not os.path.exists(fq + ".gz"):
            cmd = "bgzip {fq}".format(fq=fq)
            subprocess.check_call(cmd, shell=True)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
