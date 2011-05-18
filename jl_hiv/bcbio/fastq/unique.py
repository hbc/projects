"""Collapse fastq files into unique reads for subsequent alignment.

This helps manage complexity on highly redundant datasets.
"""
import subprocess

from bcbio.utils import (memoize_outfile, tmpfile)

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
