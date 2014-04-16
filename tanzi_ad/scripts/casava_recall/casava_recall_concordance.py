"""
Perform concordance checks between GEMINI output and a directory of known VCFs
"""
import glob
import os
import sys
import vcf
from argparse import ArgumentParser
from concordance import  concordance_calculator

def _get_vcfs(vcf_dir):
    vcfs = list(glob.glob(os.path.join(vcf_dir, "*.vcf.gz")))
    return [vcf for vcf in vcfs if "passonly" not in vcf]

def _gemini_iterator(in_file):
    with open(in_file) as in_handle:
        header = in_handle.next().split("\t")
        for line in in_handle:
            yield dict(zip(header, line.split("\t")))

def _get_vcf_handles(vcf_files):
    return [vcf.Reader(filename=v) for v in vcf_files]

def _samples_to_consider(called_families, famfile):
    """
    only consider samples from families that have been called
    """
    to_consider = []
    with open(famfile) as in_handle:
        for line in in_handle:
            family = line.split("\t")[0]
            sample_name = line.split("\t")[1]
            if family in called_families:
                to_consider.append(sample_name)
    return to_consider

if __name__ == "__main__":
    description = ("Perform concordance checks between CASAVA and Illumina "
                   "calls.")
    parser = ArgumentParser(description=description)
    parser.add_argument("gemini", help="GEMINI output. Must include a "
                "header and have the chrom, start, end, ref, alt, "
                "and variant_samples columns in the output.")
    parser.add_argument("vcf_dir", help="Directory of VCFs to scan. "
                        "Must be tabix indexed.")
    parser.add_argument("fam", help="FAM file.")
    args = parser.parse_args()

    vcf_files = _get_vcfs(args.vcf_dir)
    called_families = [os.path.basename(x).split("-")[0] for x in vcf_files]
    to_consider = _samples_to_consider(called_families, args.fam)
    vcf_handles = _get_vcf_handles(vcf_files)

    print ",".join(["chrom", "start", "end", "ref", "alt",
                    "num_concordant", "num_casava_only",
                    "num_recalled_only"])
    with open(args.gemini) as in_handle:
        header = in_handle.next().split("\t")
        for line in in_handle:
            sys.stderr.write(".")
            line = dict(zip(header, line.split("\t")))
            print concordance_calculator(line, vcf_handles, to_consider)
