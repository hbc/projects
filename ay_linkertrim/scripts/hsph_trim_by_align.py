#!/usr/bin/env python
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
from bcbio.hbc.shrna.target import identify_targets

def main(system_config_file, cur_config_file):
    config = utils.merge_config_files([system_config_file, cur_config_file])
    run_module = "bcbio.hbc.linker"
    trim_vals = config["algorithm"]["simple_trims"]
    fastq_dir = utils.add_full_path(config["dir"]["fastq"])
    cur_files = [os.path.join(fastq_dir, f) for f in config["files"]]
    dirs = {"config": utils.add_full_path(os.path.dirname(system_config_file)),
            "work" : os.getcwd(),
            "align": utils.add_full_path(config["dir"]["align"])}
    dirs["galaxy"] = os.path.dirname(utils.add_full_path(config["galaxy_config"], dirs["config"]))
    config["dir"]["trim"] = utils.add_full_path(config["dir"]["work_trim"])
    config["dir"]["fastq"] = fastq_dir
    config["dir"]["work_fastq"] = utils.add_full_path(config["dir"]["work_fastq"])
    run_parallel = parallel_runner(run_module, dirs, config, system_config_file)
    aligned = []
    for i in range(len(trim_vals.values()[0])):
        print cur_files
        in_args = [(f, i, trim_vals, config) for f in cur_files]
        align_trimmed_files = run_parallel("trim_with_aligner", in_args)
        cur_files = [x["unaligned"] for x in align_trimmed_files if x["unaligned"]]
        aligned.append([x["aligned"] for x in align_trimmed_files])
    trimmed_fastq = combine_aligned(aligned, config)
    align_bams = do_alignment(trimmed_fastq, config, dirs, run_parallel)
    count_files = count_targets(align_bams, config)

def count_targets(align_bams, config):
    """Generate count files associated with shRNA targets.
    """
    bed_targets = identify_targets(align_bams, config)

def do_alignment(trimmed_fastq, config, dirs, run_parallel):
    def _base_fname(x):
        return os.path.splitext(os.path.basename(x.replace(" ", "_")))[0]
    def _safe_fname(x):
        safe_x = x.replace(" ", "_")
        if x != safe_x and not os.path.exists(safe_x):
            os.symlink(x, safe_x)
        return safe_x
    info = {"genome_build": config["algorithm"]["genome_build"]}
    xs = run_parallel("hbc_process_alignment",
                      ((_safe_fname(f), None, info, _base_fname(f), "", dirs, config)
                       for f in trimmed_fastq))
    return [x["out_bam"] for x in xs]

def combine_aligned(aligned, config):
    """Combine aligned sequences into final output files.
    """
    trimmed = []
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
        trimmed.append(out_fname)
    return trimmed

if __name__ == "__main__":
    main(*sys.argv[1:])
