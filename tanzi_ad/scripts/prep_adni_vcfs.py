#!/usr/bin/env python
#from bcbio.variation import combine_variant_files
from argparse import ArgumentParser
import glob
import os
import yaml
from joblib import Parallel, delayed
from bcbio.variation import vcfutils, effects, population


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
    pass

def find_vcf_files(dirname):
    vcf_extensions = ("*.vcf", "*.VCF")
    queries = [os.path.join(dirname, x) for x in vcf_extensions]
    vcf_files = []
    for query in queries:
        vcf_files.extend(glob.glob(query))
    return vcf_files

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

def main(dirname, config, cores):
    vcf_files = find_vcf_files(dirname)
    prepped_files = prep_vcf_files(vcf_files, cores, config)
    merged_file = merge_vcf_files(prepped_files, cores, config)
    effects_file = effects.snpeff_effects({"vrn_file": merge_file,
                                           "genome_resources": {"aliases" : {"snpeff": "GRCh37"}},
                                           "genome_build": "GRCh37",
                                           "config": config})
    snpeff_file = effects.snpeff_effects(merged_file)
    gemini_db = population.prep_gemini_db(snpeff_file,
                                          [os.path.splitext(config["outputs"]["merge"])[0], "casava"],
                                          [{"config": config, "work_bam": "yes", "genome_build": "GRCh37"}],
                                          data)[0][1]["db"]

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
