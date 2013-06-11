from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir
import csv
import os
from bipy.utils import _download_ref, append_stem, replace_suffix, prepare_ref_file
from bipy.toolbox import (fastqc, sickle, cutadapt_tool, tophat, htseq_count,
                           rseqc, sam)
import glob
from itertools import product
from bcbio.broad import BroadRunner, picardrun
import sh


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


def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" % (stage, curr_files))


def _sam_to_bam(in_file):
    import sh
    from bipy.utils import replace_suffix
    from bcbio.utils import file_exists
    bam_file = replace_suffix(in_file, "bam")
    if file_exists(bam_file):
        return bam_file
    sh.samtools.view("-Sb", in_file, "-o", bam_file)
    return bam_file


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    stage_dict = {"download_encode": _download_encode,
                  "fastqc": _run_fastqc}

    curr_files = config["encode_file"]

    results_dir = config["dir"].get("results", "results")

    for cell_type in config["cell_types"]:
        cell_type_dir = os.path.join(results_dir, cell_type)
        safe_makedir(cell_type_dir)
        config["dir"]["results"] = cell_type_dir
        in_files = glob.glob(os.path.join(config["dir"]["data"],
                                          cell_type, "*"))
        curr_files = in_files
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

            if stage == "tophat":
                _emit_stage_message(stage, curr_files)
                tophat_config = _get_stage_config(config, stage)
                tophat_args = zip(*product(curr_files, [None], [config["ref"]],
                                           ["tophat"], [config]))
                tophat_outputs = view.map(tophat.run_with_config, *tophat_args)

                picard = BroadRunner(config["program"]["picard"])
                # convert to bam
                #args = zip(*product([picard], tophat_outputs))
                #bamfiles = view.map(picardrun.picard_formatconverter,
                #                    *args)
                bamfiles = view.map(sam.sam2bam, tophat_outputs)
                sorted_bf = view.map(sam.bamsort, bamfiles)
                view.map(sam.bamindex, sorted_bf)
                curr_files = sorted_bf

            if stage == "rseqc":
                _emit_stage_message(stage, curr_files)
                rseqc_config = _get_stage_config(config, stage)
                rseq_args = zip(*product(curr_files, [config]))
                view.map(rseqc.bam2bigwig, *rseq_args, block=False)
                view.map(rseqc.bam_stat, *rseq_args, block=False)
                view.map(rseqc.clipping_profile, *rseq_args, block=False)
                view.map(rseqc.genebody_coverage, *rseq_args, block=False)
                view.map(rseqc.junction_annotation, *rseq_args, block=False)
                view.map(rseqc.junction_saturation, *rseq_args, block=False)
                RPKM_count_files = view.map(rseqc.RPKM_count,
                                            *rseq_args)
                dirs_to_process = list(set(map(os.path.dirname,
                                               RPKM_count_files)))
                logger.info("Count files: %s" % (RPKM_count_files))
                logger.info("dirnames to process: %s" % (dirs_to_process))
                RPKM_merged = view.map(rseqc.merge_RPKM, dirs_to_process)

                view.map(rseqc.RPKM_saturation, *rseq_args, block=False)
                curr_files = tophat_outputs

            if stage == "htseq-count":
                _emit_stage_message(stage, curr_files)
                htseq_config = _get_stage_config(config, stage)
                htseq_args = zip(*product(curr_files, [config], [stage]))
                htseq_outputs = view.map(htseq_count.run_with_config,
                                         *htseq_args)
                column_names = in_files
                out_file = os.path.join(config["dir"]["results"], stage,
                                        cell_type + ".combined.counts")
                combined_out = htseq_count.combine_counts(htseq_outputs,
                                                          column_names,
                                                          out_file)
                rpkm = htseq_count.calculate_rpkm(combined_out,
                                                  config["annotation"]["file"])
                rpkm_file = os.path.join(config["dir"]["results"], stage,
                                         cell_type + ".rpkm.txt")
                rpkm.to_csv(rpkm_file, sep="\t")

            if stage == "coverage":
                logger.info("Calculating RNASeq metrics on %s." % (curr_files))
                nrun = len(curr_files)
                ref = prepare_ref_file(config["stage"][stage]["ref"],
                                              config)
                ribo = config["stage"][stage]["ribo"]
                picard = BroadRunner(config["program"]["picard"])
                out_dir = os.path.join(config["dir"]["results"], stage)
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

    # end gracefully, wait for jobs to finish, then exit
    view.wait()
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
