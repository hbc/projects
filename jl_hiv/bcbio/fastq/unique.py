"""Collapse fastq files into unique reads for subsequent alignment.

This helps manage complexity on highly redundant datasets.
"""
import subprocess

from bcbio.utils import (memoize_outfile, tmpfile)

def uniquify_reads(in_fastq, config):
    prog = config["program"]["uniquify"]
    if "bioplayground" in prog:
        return uniquify_bioplayground(in_fastq, config)
    elif "bloom" in prog:
        return uniquify_bloom(in_fastq, config)
    else:
        raise ValueError("Do not know how to uniquify with {}".format(prog))

@memoize_outfile("-unique.txt")
def uniquify_bioplayground(in_fastq, config, out_file):
    """Uniquify fastq file using Brent Pedersen's C++ fastq program.

    https://github.com/brentp/bio-playground/tree/master/reads-utils
    """
    cl = [config["program"]["uniquify"], "filter",
          "--n", str(config["algorithm"]["allowed_ns"]),
          "--unique", in_fastq]
    with open(out_file, "w") as out_handle:
        subprocess.check_call(cl, stdout=out_handle)
    return out_file

@memoize_outfile("-unique.txt")
def uniquify_bloom(fname, config, out_file):
    """Subset reads to unique set of distinct reads with bloom filter.

    https://github.com/brentp/pybloomfaster
    """
    subprocess.check_call("{prog} < {inf} > {outf}".format(prog=config["program"]["uniquify"],
                                                           inf=fname,
                                                           outf=out_file),
                          shell=True)
    return out_file
