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
from bcbio.distributed.messaging import parallel_runner

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
    run_parallel = parallel_runner(run_module, dirs, config, system_config_file)
    aligned = []
    while len(trim_vals) > 0:
        print cur_files
        five_trim, three_trim = trim_vals.pop(0)
        in_args = [(f, five_trim, three_trim, config) for f in cur_files]
        align_trimmed_files = run_parallel("trim_with_aligner", in_args)
        cur_files = [x["unaligned"] for x in align_trimmed_files]
        aligned.append([x["aligned"] for x in align_trimmed_files])
    combine_aligned(aligned, config)

def combine_aligned(aligned, config):
    """Combine aligned sequences into final output files.
    """
    out_dir = utils.safe_makedir(utils.add_full_path(config["dir"]["final"]))
    for i, fname in enumerate(config["files"]):
        # write to output file
        out_fname = os.path.join(out_dir, "{0}-trim.fastq".format(
            os.path.splitext(os.path.basename(fname))[0]))
        if not utils.file_exists(out_fname):
            with open(out_fname, "w") as out_handle:
                for in_fname in [xs[i] for xs in aligned]:
                    with open(in_fname) as in_handle:
                        out_handle.writelines(in_handle)

if __name__ == "__main__":
    main(*sys.argv[1:])
