import yaml
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import setup_logging
from rkinf.toolbox import (fastqc, cutadapt_tool, novoindex, novoalign,
                           htseq_count, blastn)
from bcbio.utils import safe_makedir
import sys
import os
from bcbio.broad import BroadRunner, picardrun
import itertools
from rkinf.utils import remove_suffix, replace_suffix
from rkinf import gtf

MAX_READ_LENGTH = 10000


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _get_short_names(input_files):
    return [remove_suffix(os.path.basename(x)) for x in
            input_files]


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
            picard = BroadRunner(config["program"]["picard"])
            args = zip(*itertools.product([picard], novoalign_outputs))
            # conver to bam
            bamfiles = view.map(picardrun.picard_formatconverter,
                                *args)
            args = zip(*itertools.product([picard], bamfiles))
            # sort bam
            sorted_bf = view.map(picardrun.picard_sort, *args)
            # index bam
            args = zip(*itertools.product([picard], sorted_bf))
            view.map(picardrun.picard_index, *args)
            curr_files = novoalign_outputs

        if stage == "htseq-count":
            nfiles = len(curr_files)
            htseq_config = _get_stage_config(config, stage)
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     curr_files,
                                     [config] * nfiles,
                                     [stage] * nfiles)
            column_names = _get_short_names(input_files)
            logger.info("Column names: %s" % (column_names))
            out_file = os.path.join(config["dir"]["results"], stage,
                                    "combined.counts")
            combined_out = htseq_count.combine_counts(htseq_outputs,
                                                      column_names,
                                                      out_file)
            rpkm = htseq_count.calculate_rpkm(combined_out,
                                              config["annotation"]["file"])
            rpkm_file = os.path.join(config["dir"]["results"], stage,
                                     "rpkm.txt")
            rpkm.to_csv(rpkm_file, sep="\t")

        if stage == "coverage":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            nrun = len(curr_files)
            ref = blastn.prepare_ref_file(config["stage"][stage]["ref"],
                                          config)
            ribo = config["stage"][stage]["ribo"]
            picard = BroadRunner(config["program"]["picard"])
            out_dir = os.path.join(results_dir, stage)
            safe_makedir(out_dir)
            out_files = [replace_suffix(os.path.basename(x),
                                        "metrics") for x in curr_files]
            out_files = [os.path.join(out_dir, x) for x in out_files]
            out_files = view.map(picardrun.picard_rnaseq_metrics,
                                 [picard] * nrun,
                                 curr_files,
                                 [ref] * nrun,
                                 [ribo] * nrun,
                                 out_files)

    stop_cluster()


if __name__ == "__main__":
    main(*sys.argv[1:])
