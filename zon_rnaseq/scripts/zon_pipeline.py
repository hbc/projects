"""
pipeline for performing RNA-seq differential expression analysis on
zebrafish
"""
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
import csv
import os
from bipy.utils import append_stem, prepare_ref_file, replace_suffix
from bipy.toolbox import (fastqc, sickle, cutadapt_tool, tophat, htseq_count,
                          deseq, annotate, rseqc, sam)
from itertools import product
import glob
from bcbio.broad import BroadRunner, picardrun


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

    input_dir = config["input_dir"]
    results_dir = config["dir"].get("results", "results")
    input_files = glob.glob(os.path.join(input_dir, "*.fq"))
    curr_files = _make_current_files(input_files)
    conditions = [os.path.basename(x).split("_")[0] for x in input_files]

    for stage in config["run"]:
        if stage == "fastqc":
            _emit_stage_message(stage, curr_files)
            fastqc_config = _get_stage_config(config, stage)
            fastqc_args = zip(*product(curr_files, [fastqc_config],
                                       [config]))
            fastqc_out = view.map(fastqc.run, *fastqc_args)
            logger.info("fastqc outfiles: %s" % (fastqc_out))

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
            # convert to bam, sort and index
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
            view.map(rseqc.RPKM_count, *rseq_args, block=False)
            view.map(rseqc.RPKM_saturation, *rseq_args, block=False)
            curr_files = tophat_outputs

        if stage == "coverage":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            nrun = len(curr_files)
            ref = prepare_ref_file(config["stage"][stage]["ref"],
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
            for comparison in deseq_config["comparisons"]:
                comparison_name = "_vs_".join(comparison)
                out_dir = os.path.join(results_dir, stage,
                                       comparison_name)
                safe_makedir(out_dir)
                indexes = [x for x, y in enumerate(conditions) if y
                           in comparison]
                htseq_files = [htseq_outputs[index] for index in indexes]
                htseq_columns = [conditions[index] for index in indexes]
                out_file = os.path.join(out_dir,
                                        comparison_name + ".counts.txt")
                combined_out = htseq_count.combine_counts(htseq_files,
                                                          htseq_columns,
                                                          out_file)
                deseq_conds = [conditions[index] for index in indexes]
                deseq_out = os.path.join(out_dir,
                                         comparison_name + ".deseq.txt")
                logger.info("Running deseq on %s with conditions %s "
                            "and writing to %s" % (combined_out,
                                                   conditions,
                                                   deseq_out))
                view.map(deseq.run, [combined_out], [deseq_conds], [deseq_out])
                annotated_file = view.map(annotate.annotate_table_with_biomart,
                                          [deseq_out],
                                          ["id"],
                                          ["ensembl_gene_id"],
                                          ["zebrafish"])


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
    from bipy.cluster import view

    main(main_config_file)
