"""
pipeline to handle performing differential analysis on ES
cells
"""
from rkinf.cluster import start_cluster, stop_cluster
import sys
import yaml
from rkinf.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
import csv
import os
from rkinf.utils import _download_ref, append_stem
from rkinf.toolbox import (fastqc, sickle, cutadapt_tool, tophat, htseq_count,
                           deseq)
from itertools import product


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


def _download_encode(input_file, config):
    """ grab the encode files they listed in their file """
    NAME_FIELD = 0
    if not os.path.exists(input_file):
        logger.info("Error %s does not exist, aborting." % (input_file))
        exit(-1)

    with open(input_file) as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        files = [x[NAME_FIELD] for x in reader]
    logger.info("Downloading %s." % (files))
    data_dir = config["dir"].get("data", "data")
    out_files = view.map(_download_ref, files, [data_dir] * len(files))

    return out_files


def _get_cell_types(input_file):
    """ keep records of the cell types for later """
    CELL_FIELD = 7
    with open(input_file) as in_handle:
        reader = csv.reader(in_handle, delimiter=";")
        cell_types = [x[CELL_FIELD] for x in reader]
    return cell_types


def _make_current_files(curr_files):
    """ makes sure the list of files is non zero and exists """
    for curr_file in curr_files:
        if not file_exists(curr_file):
            logger.error("%s does not exist or is size 0. Aborting."
                         % (curr_file))
            exit(1)
    return curr_files


def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" %(stage, curr_files))


def main(config_file):
    """ this assumes that we are keeping the same order of the files
    throughout """
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    input_dict = config["input"]
    curr_files = _make_current_files(input_dict.keys())
    input_meta = input_dict.values()

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
            curr_files = _make_current_files(cutadapt_outputs)

        if stage == "tophat":
            _emit_stage_message(stage, curr_files)
            tophat_config = _get_stage_config(config, stage)
            tophat_args = zip(*product(curr_files, [None], [config["ref"]],
                                       ["tophat"], [config]))
            tophat_outputs = view.map(tophat.run_with_config, *tophat_args)
            curr_files = _make_current_files(tophat_outputs)

        if stage == "htseq-count":
            _emit_stage_message(stage, curr_files)
            htseq_config = _get_stage_config(config, stage)
            htseq_args = zip(*product(curr_files, [config], [stage]))
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     *htseq_args)
            combined_out = os.path.join(config["dir"]["results"], stage,
                                        "all_combined.counts")
            combined_out = htseq_count.combine_counts(htseq_outputs, None,
                                                      out_file=combined_out)

        if stage == "deseq":
            _emit_stage_message(stage, curr_files)
            deseq_config = _get_stage_config(config, stage)
            out_dir = os.path.join(config["dir"]["results"], stage)
            safe_makedir(out_dir)
            for test in deseq_config["tests"]:
                indexes = [_find_file_index_for_test(input_meta,
                                                     condition) for
                                                     condition in test]
                files = [htseq_outputs[x] for x in indexes]
                conditions = [input_meta[x]["condition"] for x in indexes]
                combined_out = os.path.join(out_dir,
                                            "_".join(conditions) +
                                            "_combined.counts")
                logger.info("Combining %s to %s." % (files, combined_out))
                count_file = htseq_count.combine_counts(files, None,
                                                        out_file=combined_out)
                out_file = os.path.join(out_dir, "_".join(conditions) +
                                        "_deseq.txt")
                logger.info("Running deseq on %s with conditions %s "
                            "and writing to %s" % (count_file,
                                                   conditions,
                                                   out_file))
                view.map(deseq.run, [count_file], [conditions], [out_file])
                #deseq.run(count_file, conditions, out_file=out_file)

    # end gracefully
    stop_cluster()


def _find_file_index_for_test(in_meta, test):
    return [x[0] for x in enumerate(in_meta) if
            x[1]["condition"] == test][0]

if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    start_cluster(startup_config)
    from rkinf.cluster import view

    main(main_config_file)
