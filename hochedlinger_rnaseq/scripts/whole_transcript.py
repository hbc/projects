import yaml
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import setup_logging
from rkinf.toolbox import (fastqc, cutadapt_tool, novoindex, novoalign,
                           htseq_count)
from bcbio.utils import safe_makedir
import sys
import os

MAX_READ_LENGTH = 10000


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    from rkinf.log import logger
    start_cluster(config)

    from rkinf.cluster import view
    input_files = [os.path.join(config["dir"]["data"], x) for x in
                   config["input"]]
    results_dir = config["dir"]["results"]

    map(safe_makedir, config["dir"].values())

    curr_files = input_files

    for stage in config["run"]:
        if stage == "fastqc":
            nfiles = len(curr_files)
            logger.info("Running %s on %s" % (stage, str(curr_files)))
            fastqc_config = _get_stage_config(config, stage)
            fastqc_outputs = view.map(fastqc.run, curr_files,
                                      [fastqc_config] * nfiles,
                                      [config] * nfiles)

        if stage == "cutadapt":
            nfiles = len(curr_files)
            cutadapt_config = _get_stage_config(config, stage)
            cutadapt_outputs = view.map(cutadapt_tool.run,
                                        curr_files,
                                        [cutadapt_config] * nfiles,
                                        [config] * nfiles)
            curr_files = cutadapt_outputs

        if stage == "novoalign":
            nfiles = len(curr_files)
            novoalign_config = _get_stage_config(config, stage)
            #db = novoindex.run(config["ref"],
            #                   _get_stage_config(config, "novoindex"),
            #                   config)
            db = config["genome"]["file"]
            novoalign_outputs = view.map(novoalign.run, curr_files,
                                         [db] * nfiles,
                                         [novoalign_config] * nfiles,
                                         [config] * nfiles)
            curr_files = novoalign_outputs

        if stage == "htseq-count":
            nfiles = len(curr_files)
            htseq_config = _get_stage_config(config, stage)
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     curr_files,
                                     [config] * nfiles,
                                     [stage] * nfiles)
            combined_out = htseq_count.combine_counts(htseq_outputs,
                                                      "combined.counts")

    stop_cluster()


if __name__ == "__main__":
    main(*sys.argv[1:])
