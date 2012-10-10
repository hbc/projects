# scripts to run data analysis on tuberous sclerosis mice
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir
import csv
import os
from bipy.utils import append_stem, combine_pairs, flatten
from bipy.toolbox import (fastqc, sickle, cutadapt_tool, tophat,
                           htseq_count, deseq)
import glob
from itertools import product, repeat


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _run_trim(curr_files, config):
    logger.info("Trimming poor quality ends from %s" % (str(curr_files)))
    nfiles = len(curr_files)
    min_length = str(config["stage"]["trim"].get("min_length", 20))
    pair = str(config["stage"]["trim"].get("pair", "se"))
    platform = str(config["stage"]["trim"].get("platform", "sanger"))
    out_dir = os.path.join(config["dir"]["results"], "trimmed")
    safe_makedir(out_dir)
    out_files = [append_stem(os.path.basename(x), "trim") for
                 x in curr_files]
    out_files = [os.path.join(out_dir, x) for x in out_files]
    out_files = view.map(sickle.run, curr_files,
                         [pair] * nfiles,
                         [platform] * nfiles,
                         [min_length] * nfiles,
                         out_files)
    return out_files


def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" % (stage, curr_files))

def _find_input_files(config):
    input_dirs = config["input_dirs"]
    """ find all of the fastq files by identifier """
    identifier = config["sample_parse"]["identifier"]
    input_files = [glob.glob(os.path.join(config["dir"]["data"],
                                          input_dir,
                                          identifier))
                                          for input_dir in input_dirs]
    return list(flatten(input_files))

def _group_input_by_condition(in_files, delimiter = "_"):
    def _add_entry(d, v):
        base = os.path.basename(v)
        k = base.split(delimiter)[1]
        d[k] = d.get(k, []) + [v]
        return d

    return reduce(_add_entry, in_files, {})


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for thesis pipeline
    input_dirs = config["input_dirs"]

    results_dir = config["dir"].get("results", "results")
    input_files = _find_input_files(config)
    conditions = _group_input_by_condition(input_files)
    logger.info("Input_files: %s" % (input_files))
    logger.info("Condition groups %s" %(conditions))
    htseq_outdict = {}

    for condition, curr_files in conditions.items():
        condition_dir = os.path.join(results_dir, condition)
        safe_makedir(condition_dir)
        config["dir"]["results"] = condition_dir

        for stage in config["run"]:
            if stage == "fastqc":
                _emit_stage_message(stage, curr_files)
                fastqc_config = _get_stage_config(config, stage)
                fastqc_args = zip(*product(curr_files, [fastqc_config],
                                           [config]))
                view.map(fastqc.run, *fastqc_args)

            if stage == "cutadapt":
                _emit_stage_message(stage, curr_files)
                cutadapt_config = _get_stage_config(config, stage)
                cutadapt_args = zip(*product(curr_files, [cutadapt_config],
                                             [config]))
                cutadapt_outputs = view.map(cutadapt_tool.run, *cutadapt_args)
                curr_files = cutadapt_outputs

            if stage == "sickle":
                _emit_stage_message(stage, curr_files)
                pairs = combine_pairs(curr_files)
                first = [x[0] for x in pairs]
                second = [x[1] for x in pairs]
                fixed = view.map(sickle.run_with_config,
                                 first, second, [config] * len(first))
                curr_files = list(flatten(fixed))


            if stage == "tophat":
                _emit_stage_message(stage, curr_files)
                tophat_config = _get_stage_config(config, stage)
                pairs = combine_pairs(curr_files)
                first = [x[0] for x in pairs]
                second = [x[1] for x in pairs]
                logger.info("first %s" % (first))
                logger.info("second %s" % (second))

                #tophat_args = zip(*product(first, second, [config["ref"]],
                #                           ["tophat"], [config]))
                tophat_outputs = view.map(tophat.run_with_config,
                                          first, second,
                                          [config["ref"]] * len(first),
                                          ["tophat"] * len(first),
                                          [config] * len(first))
                curr_files = tophat_outputs

            if stage == "htseq-count":
                _emit_stage_message(stage, curr_files)
                htseq_config = _get_stage_config(config, stage)
                htseq_args = zip(*product(curr_files, [config], [stage]))
                htseq_outputs = view.map(htseq_count.run_with_config,
                                         *htseq_args)
                htseq_outdict[condition] = htseq_outputs

    # combine htseq-count files and run deseq on them
    conditions = []
    htseq_files = []
    for condition, files in htseq_outdict.items():
        conditions += list(repeat(condition, len(files)))
        htseq_files += files

    _emit_stage_message("deseq", htseq_files)
    deseq_dir = os.path.join(results_dir, "deseq")
    safe_makedir(deseq_dir)
    out_file = os.path.join(deseq_dir, "all.combined.counts")
    combined_out = htseq_count.combine_counts(htseq_outputs,
                                              files,
                                              out_file)
    out_file = os.path.join(deseq_dir, "deseq.txt")
    logger.info("Running deseq on %s with conditions %s and "
                " writing to %s." % (combined_out, conditions, out_file))
    deseq.run(combined_out, conditions, out_file=out_file)

    # end gracefully
    stop_cluster()

if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    start_cluster(startup_config)
    from bipy.cluster import view

    main(main_config_file)
