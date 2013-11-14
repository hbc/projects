#!/usr/bin/env python
from argparse import ArgumentParser
import glob
import os
import yaml
from joblib import Parallel, delayed
import subprocess
from itertools import dropwhile
from six import advance_iterator
import sys


from bcbio.variation import vcfutils, effects, population
from bcbio.distributed.messaging import parallel_runner
from bcbio.distributed.transaction import file_transaction


# def get_sample_name(vcf_file):
#     with open(vcf_file) as in_handle:
#         in_handle = dropwhile(lambda x: not x.startswith("#CHROM"), in_handle)
#         l = advance_iterator(in_handle)
#         name = l.split()[9].strip()
#         return name.split("_")[0]

def get_sample_name(vcf_file):
    base, _ = os.path.splitext(vcf_file)
    return os.path.basename(base.split("_SNPs")[0].upper())

def prep_vcf_files(vcf_files, cores, config):
    out_files = Parallel(cores)(delayed(run_illumina_prep)(s, config) for s in vcf_files)
    return out_files


def merge_vcf_files(sample_files, cores, config):
    out_file = config["outputs"]["merge"]
    config["algorithm"] = {}
    run_parallel = parallel_runner({"type": "local", "cores": min(cores, 8)}, {}, config)
    vcfutils.parallel_combine_variants(sample_files, out_file, config["ref"]["GRCh37"],
                                       config, run_parallel)
    return out_file


def find_vcf_files(dirname):
    vcf_extensions = ("*.vcf", "*.VCF")
    queries = [os.path.join(dirname, x) for x in vcf_extensions]
    vcf_files = []
    for query in queries:
        vcf_files.extend(glob.glob(query))
    return vcf_files

def run_illumina_prep(vcf_file, config):
    tmp_dir = config.get("tmpdir", os.getcwd())
    if not os.path.exists(tmp_dir):
        try:
            os.makedirs(tmp_dir)
        except OSError:
            assert os.path.exists(tmp_dir)
    sample_id = get_sample_name(vcf_file)
    out_file = os.path.join(os.getcwd(), "%s.vcf" % sample_id)
    if not os.path.exists(out_file):
        print sample_id, vcf_file, out_file
        subprocess.check_call(["java", "-Xms1g", "-Xmx2g", "-jar", config["bcbio.variation"],
                               "variant-utils", "illumina", vcf_file,
                               sample_id, config["ref"]["GRCh37"],
                               config["ref"]["hg19"],
                               "--outdir", os.getcwd(),
                               "--tmpdir", tmp_dir])
    return out_file

def load_gemini_db(vcf_file, ped_file, cores):
    cmd = ("gemini load --passonly -v {vcf_file} -p {ped_file} -t snpEff --cores "
           "{cores} {tx_gemini_db}")
    base, _ = os.path.splitext(vcf_file)
    gemini_db = base + ".db"
    with file_transaction(gemini_db) as tx_gemini_db:
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return gemini_db

def main(dirname, config, cores):
    vcf_files = find_vcf_files(dirname)
    prepped_files = prep_vcf_files(vcf_files, cores, config)
    merged_file = merge_vcf_files(prepped_files, cores, config)
    effects_file = effects.snpeff_effects({"vrn_file": merged_file,
                                           "genome_resources": {"aliases" : {"snpeff": "GRCh37"}},
                                           "genome_build": "GRCh37",
                                           "config": config})

    gemini_db = load_gemini_db(effects_file, config["ped"], cores)


#    gemini_db = population.prep_gemini_db(effects_file, call_id, samples, data)[0][1]["db"]


if __name__ == "__main__":
    parser = ArgumentParser(description="Clean, merge and load a directory of "
                            "VCF files into GEMINI.")
    parser.add_argument("config", help="Configuration YAML file.")
    parser.add_argument("vcf_dir", help="Directory of vcf files to process.")
    parser.add_argument("--cores", default=1, type=int, help="Cores to use.")
    args = parser.parse_args()
    cores = 1
    with open(args.config) as config_handle:
        config = yaml.load(config_handle)
    main(args.vcf_dir, config, args.cores)
