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
        sample_name = os.path.basename(exp["files"][0].replace("Jon_Li_", "S")).split("_")[0]
    base_out = os.path.join(out_dir, "%s-%s" % (exp_name, sample_name))

    final = "%s.fq.gz" % base_out
    fq1 = "%s-1.fq" % base_out
    fq2 = "%s-2.fq" % base_out
    if not os.path.exists(final):
        ref_file = "%s.ndx" % os.path.splitext(exp["ref"])[0]
        fastq_file, pair_file = exp["files"]
        for orig, out in zip(exp["files"], [fq1, fq2]):
            cmd = ("novoalign -c 32 -o SAM -d {ref_file} -f {orig} "
                   "| samtools view -b -S -u -F 4 - "
                   "| bamtofastq F=/dev/null F2=/dev/null S={out} O=/dev/null O2=/dev/null")
            if not os.path.exists(out):
                subprocess.check_call(cmd.format(**locals()), shell=True)
        cmd = "cat {fq1} {fq2} | bgzip -c > {final}"
        subprocess.check_call(cmd.format(**locals()), shell=True)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
