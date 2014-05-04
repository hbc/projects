#!/usr/bin/env python
"""Run squaring off batch scripts for priority 1 and 2 families.
"""
import csv
import glob
import os
import shutil
import subprocess
import sys

import joblib

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.variation import vcfutils

priority_file = "/n/hsphS10/hsphfs1/chb/projects/tanzi_ad/data/gwas/AD-Master-v2.csv"
input_vcf_dir = "/n/hsphS10/hsphfs1/chb/projects/tanzi_ad/data/recall_variants/freebayes"
bam_dir = "/n/hsphS10/hsphfs2/tanzi_recalled"
name = "tanzi_ad_p1and2-square"
ref_file = "/n/hsphS10/hsphfs1/chb/biodata/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"

def main(cores=1):
    start_dir = os.getcwd()
    work_dir = utils.safe_makedir("/scratch/square")
    priorities = set(["1", "2"])
    list_file = get_input_list(start_dir, priorities)
    ensure_bam_index(list_file)
    # Ensure input CRAMs are indexed; gets IO bound quickly so limit cores
    cram_cores = min(int(cores), 6)
    for cindex in joblib.Parallel(cram_cores)(joblib.delayed(index_cram)(x) for x in find_crams(list_file)):
        print cindex
    with utils.chdir(work_dir):
        out_file = run_squaring(list_file, name, ref_file, cores)
    for ext in ["", ".tbi"]:
        new_file = os.path.join(start_dir, os.path.basename(out_file) + ext)
        if not utils.file_exists(new_file):
            shutil.copy(out_file + ext, new_file)

def run_squaring(list_file, name, ref_file, cores):
    mem = 4 * int(cores)
    mem_opts = ["-Xms%sG" % (mem // 2), "-Xmx%sG" % mem]
    out_file = os.path.join(os.getcwd(), "%s.vcf.gz" % name)
    if not utils.file_exists(out_file):
        subprocess.check_call(["bcbio-variation-recall"] + mem_opts +
                              ["square", out_file, ref_file, list_file,
                               "--caller", "freebayes", "--cores", str(cores)])
    return out_file

def ensure_bam_index(in_file):
    with open(in_file) as in_handle:
        for line in (l.strip() for l in in_handle):
            if line.endswith(".bam"):
                if not os.path.exists(line + ".bai"):
                    old_index = "%s.bai" % os.path.splitext(line)[0]
                    if os.path.exists(old_index):
                        shutil.move(old_index, line + ".bai")
                assert os.path.exists(line + ".bai"), line

def index_cram(fname):
    out_file = "%s.crai" % fname
    if not utils.file_exists(out_file):
        print "Indexing", fname
        with file_transaction(out_file) as tx_out_file:
            tx_in_file = os.path.splitext(tx_out_file)[0]
            utils.symlink_plus(fname, tx_in_file)
            subprocess.check_call(["cram_index", tx_in_file])
    return out_file

def find_crams(in_file):
    with open(in_file) as in_handle:
        for line in (l.strip() for l in in_handle):
            if line.endswith(".cram"):
                yield line

def get_input_list(work_dir, priorities):
    list_file = os.path.join(work_dir, "p1and2-input-files.txt")
    if not utils.file_exists(list_file):
        families = read_families_by_priority(priority_file, priorities)
        vcf_files = [get_vcf(input_vcf_dir, fam) for fam in families]
        bam_files = []
        for vcf_file in vcf_files:
            bam_files.extend(get_bams(vcf_file, bam_dir))

        with open(list_file, "w") as out_handle:
            for fname in vcf_files + bam_files:
                out_handle.write(fname + "\n")
    return list_file

def get_vcf(vcf_dir, fam):
    vcfs = sorted(glob.glob(os.path.join(vcf_dir, "%s-*vcf*" % fam)))
    for ending in [".vcf.gz", ".vcf"]:
        for f in vcfs:
            if f.endswith(ending):
                return f
    raise ValueError("Did not find VCF for %s in %s" % (fam, vcf_dir))

def get_bams(vcf_file, bam_dir):
    out = []
    for sample in vcfutils.get_samples(vcf_file):
        bam_files = glob.glob(os.path.join(bam_dir, "*", "final", sample, "%s-*am" % sample))
        assert len(bam_files) > 0, "Did not find BAM files for %s: %s" % (sample, bam_files)
        if len(bam_files) > 1:
            bam_files = [x for x in bam_files if x.endswith(".bam")]
        out.append(bam_files[0])
    return out

def read_families_by_priority(fname, priorities):
    families = set([])
    with open(fname) as in_handle:
        reader = csv.reader(in_handle)
        reader.next()  # header
        for parts in reader:
            _, family_id, _, priority = parts[:4]
            status_flag = parts[16]
            if status_flag != "Exclude" and priority.strip() in priorities:
                families.add(family_id)
    return sorted(list(families))


if __name__ == "__main__":
    main(*sys.argv[1:])
