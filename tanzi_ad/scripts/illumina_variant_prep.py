#!/usr/bin/env python
"""Prepare Illumina called variants, merging into single sample GATK-compatible VCFs.

Usage:
  illumina_variant_prep.py <in config> <cores>
"""
import csv
import glob
import pprint
import os
import shutil
import subprocess
import sys

import yaml

import joblib

from bcbio import broad, utils
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
    config["algorithm"]["num_cores"] = cores
    data = {"config": config, "dirs": {"work": os.getcwd()}, "name": [""]}
    prepare_plink_vcftools(noexclude_file, config)
    gemini_db = population.prep_gemini_db([noexclude_file],
                                          [os.path.splitext(config["outputs"]["merge"])[0], "casava"],
                                          [{"config": config, "work_bam": "yes", "genome_build": "GRCh37"}],
                                          data)[0][1]["db"]
    print gemini_db

def merge_vcf_files(sample_files, cores, config):
    out_file = config["outputs"]["merge"]
    config["algorithm"] = {}
    run_parallel = parallel_runner({"type": "local", "cores": min(cores, 8)}, {}, config)
    vcfutils.parallel_combine_variants(sample_files, out_file, config["ref"]["GRCh37"],
                                       config, run_parallel)
    return out_file

def _remove_plink_problems(in_vcf):
    """Remove lines which cause issues feeding into plink.
    """
    chr_remap = {"X": "23", "Y": "24"}
    out_vcf = "%s-plinkready%s" % os.path.splitext(in_vcf)
    support_chrs = set(["G", "A", "T", "C", ","])
    if not utils.file_exists(out_vcf):
        with open(in_vcf) as in_handle:
            with open(out_vcf, "w") as out_handle:
                for line in in_handle:
                    # Mitochondrial calls are not called with
                    # same reference, causing plink conversion errors
                    if line.startswith("MT"):
                        line = None
                    # Remap sex chromosomes to plink binary numbering
                    elif line.startswith(tuple(chr_remap.keys())):
                        line = chr_remap[line[0]] + line[1:]
                    elif not line.startswith("#"):
                        parts = line.split("\t", 7)
                        if parts[6] not in set(["PASS", "."]):
                            line = None
                        else:
                            for to_check in parts[3:5]:
                                # Exclude non-GATC inputs
                                if len(set(to_check) - support_chrs) > 0:
                                    line = None
                                    break

                    if line:
                        out_handle.write(line)
    return out_vcf

def _bgzip_vcf(in_vcf):
    out_file = "%s.gz" % in_vcf
    if not utils.file_exists(out_file):
        subprocess.check_call(["bgzip", in_vcf])
        #with open(in_vcf, "w") as out_handle:
        #    out_handle.write("Gzipped to %s" % out_file)
    return out_file

def vcftools_vcf_to_tped(in_vcf, base_dir, chromosome, config):
    """Use VCFtools to convert a single chromosome into a tped file.
    """
    out_base = os.path.join(base_dir, "%s-%s" % (config["outputs"]["plink"], chromosome))
    out_tped = out_base + ".tped"
    if not utils.file_exists(out_tped):
        subprocess.check_call(["vcftools", "--plink-tped", "--gzvcf", in_vcf,
                               "--chr", str(chromosome), "--out", out_base])
    return out_tped

def _merge_tped_tfam(fnames, out_base, config):
    """Merge all tped files together and copy over fam lines by samples to tfam.
    """
    tped_file = "%s-transposed.tped" % out_base
    if not utils.file_exists(tped_file):
        with open(tped_file, "w") as out_handle:
            for fname in fnames:
                with open(fname) as in_handle:
                    for line in in_handle:
                        out_handle.write(line)
    tfam_file = "%s.tfam" % os.path.splitext(tped_file)[0]
    if not utils.file_exists(tfam_file):
        fam_lines = {}
        with open(config["array"]["fam"]) as in_handle:
            for line in in_handle:
                fam_lines[line.split()[1]] = line
        with open(tfam_file, "w") as out_handle:
            base_tfam = "%s.tfam" % os.path.splitext(fnames[0])[0]
            with open(base_tfam) as in_handle:
                for line in in_handle:
                    out_line = fam_lines[line.split()[1]]
                    out_handle.write(out_line)
    return tped_file, tfam_file

def _plink_transposed_to_bed(tped_file, tfam_file, out_base):
    """Generate finalize binary PED files from tped and tfam inputs.
    """
    out_file = "%s.bed"
    if not utils.file_exists(out_file):
        subprocess.check_call(["plink", "--tped", tped_file, "--tfam", tfam_file,
                               "--make-bed", "--out", out_base])
    return out_file

def prepare_plink_vcftools(in_vcf, config):
    """Prepare binary PED files using vcftools and plink.
    """
    clean_vcf = _remove_plink_problems(in_vcf)
    bgzip_vcf = _bgzip_vcf(clean_vcf)
    out_base = os.path.join(os.path.dirname(in_vcf), config["outputs"]["plink"])
    parts_dir = utils.safe_makedir(out_base + "-chromosomes")
    cores = config["algorithm"]["num_cores"]
    out_files = list(joblib.Parallel(cores)(joblib.delayed(vcftools_vcf_to_tped)
                                            (bgzip_vcf, parts_dir, c, config) for c in range(1, 25)))
    tped_file, tfam_file = _merge_tped_tfam(out_files, out_base, config)
    return _plink_transposed_to_bed(tped_file, tfam_file, out_base)

def prepare_plink_gatk(in_vcf, config):
    """Prepare binary PED files for input to plink using GATK walkers.

    XXX Does not produce valid binary PED files for large numbers of input samples.
    """
    clean_vcf = _remove_plink_problems(in_vcf)
    runner = broad.runner_from_config(config)
    out_base = os.path.join(os.path.dirname(in_vcf), config["outputs"]["plink"])
    args = ["-T", "VariantsToBinaryPed", "-R", config["ref"]["GRCh37"],
            "--variant", clean_vcf, "--minGenotypeQuality", "0",
            "--metaData", config["array"]["fam"],
            "--bed", out_base + ".bed", "--bim", out_base + ".bim", "--fam", out_base + ".fam"]
    runner.run_gatk(args)
    print out_base

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
