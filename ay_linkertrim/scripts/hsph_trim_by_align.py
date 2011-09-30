#!/usr/bin/env
"""Trim linker sequences using alignments, prioritizing trimmings that map.

This is useful for cases where the linker sizes are variable. Without
a defined pattern to fall back on, try variable trimming options and maintain
mapping reads.

Usage:
  hsph_trim_by_align.py <system YAML config> <YAML config>
"""
import os
import sys

import yaml

from bcbio import utils
from bcbio.distributed.messaging import run_parallel

def main(system_config_file, cur_config_file):
    config = utils.merge_config_files([system_config_file, cur_config_file])
    run_module = "bcbio.hbc.linker"
    trim_vals = config["algorithm"]["simple_trims"]
    fastq_dir = utils.add_full_path(config["dir"]["fastq"])
    cur_files = [os.path.join(fastq_dir, f) for f in config["files"]]
    dirs = {"config": utils.add_full_path(os.path.dirname(system_config_file)),
            "work" : os.getcwd()}
    config["dir"]["trim"] = utils.add_full_path(config["dir"]["trim"])
    config["dir"]["fastq"] = fastq_dir
    aligned = []
    while len(trim_vals) > 0:
        print cur_files
        five_trim, three_trim = trim_vals.pop(0)
        in_args = [(f, five_trim, three_trim, config) for f in cur_files]
        align_trimmed_files = run_parallel(run_module, "trim_with_aligner",
                                           in_args, dirs,
                                           config, system_config_file)
        print align_trimmed_files
        cur_files = [x["unaligned"] for x in align_trimmed_files]
        aligned.append([x["aligned"] for x in align_trimmed_files])
    print aligned

if __name__ == "__main__":
    main(*sys.argv[1:])
