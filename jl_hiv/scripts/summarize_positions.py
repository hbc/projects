#!/usr/bin/env python
"""Summarize constant and variable positions of interest in tab delimited files.

This prepares positions of interest for detailed interrogation with Hadoop
scripts in snp-assess. The output is two tab delimited files: location_constant.tsv
and location_variable.tsv that can be fed into examinations of off-target and minimum
coverage questions.

Usage:
  summarize_positions.py <YAML config>
"""
import os
import csv
import sys
import operator
from contextlib import nested

import yaml
from Bio import SeqIO

from bcbio.varcall.mixed import read_call_file

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    chr_name = get_chr_name(config["ref"])
    constant_file, variable_file = get_out_files(config)
    with nested(open(constant_file, "w"), open(variable_file, "w")) as \
         (constant_handle, variable_handle):
        constant_writer = csv.writer(constant_handle, dialect="excel-tab")
        variable_writer = csv.writer(variable_handle, dialect="excel-tab")
        for e in config["expected"]:
            constant, variable = constant_variable_positions(e["file"], e["offset"])
            write_out(constant_writer, constant, chr_name)
            write_out(variable_writer, variable, chr_name)

def constant_variable_positions(in_file, offset):
    constant = []
    variable = []
    for (_, pos), bases in read_call_file(in_file, ignore_space=True):
        exp_pos = pos + offset
        exp_bases = [(b, f) for b, f in bases.iteritems() if f is not None]
        if len(exp_bases) == 1:
            cur_base, cur_freq = exp_bases[0]
            constant.append((exp_pos, cur_base, cur_freq))
        else:
            exp_bases.sort(key=operator.itemgetter(1))
            cur_base, cur_freq = exp_bases[0]
            variable.append((exp_pos, cur_base, cur_freq))
    return constant, variable

def write_out(writer, info, chr_name):
    for pos, base, freq in info:
        writer.writerow([chr_name, pos, base, freq])

def get_out_files(config):
    out_dir = os.path.join(os.getcwd(), config["dir"]["vrn"], "hadoop")
    files = []
    for name in ["constant", "variable"]:
        cur_dir = os.path.join(out_dir, name)
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)
        files.append(os.path.join(cur_dir, "location_{0}.tsv".format(name)))
    return files

def get_chr_name(in_file):
    name = None
    with open(in_file) as in_handle:
        for rec in SeqIO.parse(in_handle, "fasta"):
            assert name is None
            name = rec.id
    assert name is not None
    return name

if __name__ == "__main__":
    main(*sys.argv[1:])
