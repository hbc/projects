"""Collapse fastq files into unique reads for subsequent alignment.

This helps manage complexity on highly redundant datasets.
"""
import subprocess
import collections

import yaml
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.utils import (memoize_outfile, tmpfile)

# Create unique reads from input FASTQ files

def uniquify_reads(in_fastq, config):
    prog = config["program"]["uniquify"]
    if "bioplayground" in prog:
        unique_file = uniquify_bioplayground(in_fastq, config)
    elif "bloom" in prog:
        unique_file = uniquify_bloom(in_fastq, config)
    else:
        raise ValueError("Do not know how to uniquify with {}".format(prog))
    min_counts = int(config["algorithm"].get("min_unique_counts", 1))
    if min_counts > 1:
        unique_file = filter_by_mincounts(unique_file, in_fastq, min_counts)
    count_file = generate_unique_counts(in_fastq, unique_file)
    return unique_file, count_file

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

@memoize_outfile("-filter.txt")
def filter_by_mincounts(unique, orig, min_counts, out_file):
    """Filter reads by minimum counts.
    """
    counts = collections.defaultdict(int)
    with open(orig) as in_handle:
        for _, seq, _ in FastqGeneralIterator(in_handle):
            counts[seq] += 1
    with open(unique) as in_handle:
        with open(out_file, "w") as out_handle:
            i = 0
            for _, seq, qual in FastqGeneralIterator(in_handle):
                if counts[seq] >= min_counts:
                    out_handle.write("@%s\n%s\n+\n%s\n" % (i, seq, qual))
                    i += 1
    return out_file

# Statistics for generated unique reads

def _blank_read_stats(name, seq):
    return {"seq": seq,
            "count": 0,
            "qualities": collections.defaultdict(lambda: collections.defaultdict(int))}

def _finalize_read_stats(stats):
    out_quals = {}
    for key, val in stats["qualities"].iteritems():
        out_quals[key] = dict(val)
    stats["qualities"] = out_quals
    return stats

@memoize_outfile("-unique-counts.yaml")
def generate_unique_counts(orig, unique, out_file):
    """Generate file of counts and quality scores for uniquified reads.
    """
    stats = {}
    print "Generating unique counts, preparing unique"
    with open(unique) as in_handle:
        for name, seq, _ in FastqGeneralIterator(in_handle):
            stats[seq] = _blank_read_stats(name, seq)
    print "Generating unique counts, reading full file"
    with open(orig) as in_handle:
        for name, seq, qual in FastqGeneralIterator(in_handle):
            if stats.has_key(seq):
                stats[seq]["count"] += 1
                #for i, q in enumerate(qual):
                #    stats[seq]["qualities"][i][q] += 1
    final = map(_finalize_read_stats, stats.itervalues())
    with open(out_file, "w") as out_handle:
        yaml.dump(final, out_handle, allow_unicode=False)
