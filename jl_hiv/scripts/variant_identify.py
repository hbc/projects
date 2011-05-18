#!/usr/bin/env python
"""Identify variants in barcoded multiplexed HIV samples.

Usage:
  variant_identify.py <YAML config>
"""
import sys

import yaml

from bcbio.utils import create_dirs
from bcbio.fastq.barcode import demultiplex, convert_illumina_oldstyle
from bcbio.fastq.unique import uniquify_bioplayground

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    print config
    create_dirs(config)
    for cur in config["input"]:
        in_fastq = cur["fastq"]
        if cur.get("old_style_barcodes", False):
            in_fastq = convert_illumina_oldstyle(in_fastq)
        bc_files = demultiplex(in_fastq, cur["barcodes"],
                               config["dir"]["tmp"], config)
        for bc_file in bc_files:
            uniquify_fastq(bc_file, config)

def uniquify_fastq(in_file, config):
    """Convert fastq file into collapsed fastq with unique reads.
    """
    unique_file = uniquify_bioplayground(in_file, config)
    print unique_file

if __name__ == "__main__":
    main(*sys.argv[1:])
