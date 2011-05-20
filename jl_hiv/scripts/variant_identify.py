#!/usr/bin/env python
"""Identify variants in barcoded multiplexed HIV samples.

Usage:
  variant_identify.py <YAML config>
"""
import os
import sys
import glob
import subprocess

import yaml
import khmer
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.utils import create_dirs, map_wrap, cpmap
from bcbio.fastq.barcode import demultiplex, convert_illumina_oldstyle
from bcbio.fastq.unique import uniquify_bioplayground
from bcbio.ngsalign import novoalign

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    print config
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
    unique_file = uniquify_bioplayground(in_file, config)
    align_sam = novoalign.align(config["dir"]["align"], ref_index, unique_file,
                                qual_format=cur_config.get("format", None))
    print align_sam
    align_bam = to_bamsort(align_sam, in_file, config, config_file)
    print align_bam
    ktable = make_ktable(in_file, config["kmer_size"])

def make_ktable(in_fastq, kmer_size):
    ktable = khmer.new_ktable(kmer_size)
    with open(in_fastq) as in_handle:
        for (_, seq, _) in FastqGeneralIterator(in_handle):
            if seq.find("N") == -1:
                ktable.consume(seq)
    return ktable

def to_bamsort(sam_file, fastq_file, config, config_file):
    sample_name = os.path.splitext(os.path.basename(sam_file))[0]
    cl = [config["program"]["bamsort"], "--name=%s" % sample_name,
          config_file, sam_file, config["ref"], fastq_file]
    subprocess.check_call(cl)
    return glob.glob("%s*-sort.bam" % os.path.splitext(sam_file)[0])

if __name__ == "__main__":
    main(*sys.argv[1:])
