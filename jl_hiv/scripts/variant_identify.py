#!/usr/bin/env python
"""Identify variants in barcoded multiplexed HIV samples.

Usage:
  variant_identify.py <YAML config>
"""
import os
import sys
import csv
import glob
import operator
import subprocess
import collections
from contextlib import closing

import yaml
import pysam
import khmer
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.utils import create_dirs, map_wrap, cpmap
from bcbio.broad import BroadRunner
from bcbio.fastq.barcode import demultiplex, convert_illumina_oldstyle
from bcbio.fastq.unique import uniquify_bioplayground
from bcbio.ngsalign import novoalign

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    ref_index = novoalign.refindex(config["ref"], kmer_size=13, step_size=1)
    create_dirs(config)
    for cur in config["input"]:
        in_fastq = cur["fastq"]
        if cur.get("old_style_barcodes", False):
            in_fastq = convert_illumina_oldstyle(in_fastq)
        bc_files = demultiplex(in_fastq, cur["barcodes"],
                               config["dir"]["tmp"], config)
        with cpmap(config["algorithm"]["cores"]) as cur_map:
            for _ in cur_map(process_fastq, ((bc_file, ref_index, cur, config, config_file)
                                             for bc_file in bc_files)):
                pass

@map_wrap
def process_fastq(in_file, ref_index, cur_config, config, config_file):
    picard = BroadRunner(config["program"]["picard"], config["program"]["gatk"],
                         config["algorithm"]["java_memory"])
    unique_file = uniquify_bioplayground(in_file, config)
    align_sam = novoalign.align(config["dir"]["align"], ref_index, unique_file,
                                qual_format=cur_config.get("format", None))
    align_bam = to_bamsort(align_sam, in_file, config, config_file)
    realign_bam = picard.run_fn("gatk_realigner", align_bam, config["ref"],
                                deep_coverage=True)

    picard.run_fn("picard_index", realign_bam)
    print realign_bam
    position_percent_file(align_bam, in_file, config)

def position_percent_file(align_bam, read_file, config):
    kmer_size = config["algorithm"]["kmer_size"]
    min_thresh = config["algorithm"]["detection_thresh"]
    bases = ["A", "C", "G", "T"]
    out_file = os.path.join(config["dir"]["vrn"], "%s-variations.tsv" %
                            os.path.splitext(os.path.basename(align_bam))[0])
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(["space", "pos"] + bases)
            ktable = make_ktable(read_file, kmer_size)
            for chrom, pos, kmers in positional_kmers(align_bam, kmer_size):
                base_percents = {}
                for base, percent in base_kmer_percents(kmers, ktable, min_thresh):
                    base_percents[base] = "%.1f" % (percent * 100.0)
                writer.writerow([chrom, pos] + [base_percents.get(b, "") for b in bases])

def base_kmer_percents(kmers, ktable, min_thresh):
    """Retrieve percentages of each base call based on k-mer counts.
    """
    kmer_counts = []
    for kmer, base in kmers:
        kmer_counts.append(ktable.get(kmer))
    total = float(sum(kmer_counts))
    base_counts = collections.defaultdict(int)
    for (kmer, base), kcount in zip(kmers, kmer_counts):
        if kcount / total > min_thresh:
            base_counts[base] += kcount
    pass_total = float(sum(base_counts.values()))
    final = []
    for base, count in base_counts.iteritems():
        final.append((base, count / pass_total))
    final.sort(key=operator.itemgetter(1), reverse=True)
    return final

def positional_kmers(in_bam, kmer_size):
    """Retrieve informative kmers at each piled up position in an alignment.
    """
    with closing(pysam.Samfile(in_bam, 'rb')) as work_bam:
        for col in work_bam.pileup():
            space = work_bam.getrname(col.tid)
            kmers = list(set(filter(lambda x: x is not None,
                                    [_read_surround_region(r, kmer_size) for r in col.pileups])))
            yield space, col.pos, kmers

def _read_surround_region(read, kmer_size):
    """Provide context for an aligned read at a particular position.

    Requires a full length kmer at each side so excludes information near
    start and end of reads.
    """
    assert kmer_size % 2 == 1, "Need odd kmer size"
    extend = (kmer_size - 1) // 2
    if read.indel == 0:
        seq = read.alignment.seq
        if read.qpos >= extend and read.qpos < len(seq) - extend:
            kmer = seq[read.qpos-extend:read.qpos+extend+1]
            assert len(kmer) == kmer_size, (kmer, seq, read.qpos, len(seq))
            call = seq[read.qpos]
            return (kmer, call)

def make_ktable(in_fastq, kmer_size):
    ktable = khmer.new_ktable(kmer_size)
    with open(in_fastq) as in_handle:
        #i = 0
        for (_, seq, _) in FastqGeneralIterator(in_handle):
            #i += 1
            #if i > 1e5: break
            if seq.find("N") == -1:
                ktable.consume(seq)
    return ktable

def to_bamsort(sam_file, fastq_file, config, config_file):
    sample_name = os.path.splitext(os.path.basename(sam_file))[0]
    cl = [config["program"]["bamsort"], "--name=%s" % sample_name,
          config_file, sam_file, config["ref"], fastq_file]
    subprocess.check_call(cl)
    return glob.glob("%s*-sort.bam" % os.path.splitext(sam_file)[0])[0]

if __name__ == "__main__":
    main(*sys.argv[1:])
