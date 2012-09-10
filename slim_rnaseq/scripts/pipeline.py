from rkinf.cluster import start_cluster, stop_cluster
import sys
import yaml
from rkinf.log import setup_logging, logger
from bcbio.utils import safe_makedir
import csv
import os
from rkinf.utils import _download_ref, append_stem
from rkinf.toolbox import fastqc, sickle, cutadapt_tool, tophat, htseq_count


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _run_fastqc(curr_files, config):
    logger.info("Running fastqc on %s" % (str(curr_files)))
    nfiles = len(curr_files)
    fastqc_config = config["stage"]["fastqc"]
    out_files = view.map(fastqc.run, curr_files,
                         [fastqc_config] * nfiles,
                         [config] * nfiles)
    return out_files


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


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    stage_dict = {"download_encode": _download_encode,
                  "fastqc": _run_fastqc}

    curr_files = config["encode_file"]

    for stage in config["run"]:
        if stage == "download_encode":
            curr_files = _download_encode(config["encode_file"], config)
        if stage == "fastqc":
            _run_fastqc(curr_files, config)
        if stage == "trim":
            _run_trim(curr_files, config)
        if stage == "cutadapt":
            nfiles = len(curr_files)
            logger.info("Running %s on %s" % (stage, str(curr_files)))
            cutadapt_config = _get_stage_config(config, stage)
            cutadapt_outputs = view.map(cutadapt_tool.run,
                                        curr_files,
                                        [cutadapt_config] * nfiles,
                                        [config] * nfiles)
            curr_files = cutadapt_outputs
        if stage == "tophat":
            nfiles = len(curr_files)
            tophat_config = _get_stage_config(config, stage)
            tophat_outputs = view.map(tophat.run_with_config, curr_files,
                                      [None] * nfiles,
                                      [config["ref"]] * nfiles,
                                      ["tophat"] * nfiles,
                                      [config] * nfiles)
            curr_files = tophat_outputs

        if stage == "htseq-count":
            nfiles = len(curr_files)
            htseq_config = _get_stage_config(config, stage)
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     curr_files
                                     [config] * nfiles,
                                     [stage] * nfiles)
            # combine the counts for each type

            combined_out = htseq_count.combine_counts(htseq_outputs, None,
                                                      "combined.counts")

    cell_types = _get_cell_types(config["encode_file"])
    logger.info("files: %s" % (curr_files))
    logger.info("types: %s" % (cell_types))

    # end gracefully
    stop_cluster()

if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    start_cluster(startup_config)
    from rkinf.cluster import view

    main(main_config_file)
