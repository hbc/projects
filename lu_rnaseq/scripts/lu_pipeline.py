# scripts to run data analysis on tuberous sclerosis mice
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
import csv
import os
from bipy.utils import (append_stem, combine_pairs, flatten, dict_to_vectors,
                        prepare_ref_file, replace_suffix)
from bipy.toolbox import (fastqc, sickle, cutadapt_tool, tophat,
                           htseq_count, deseq, fastq, annotate, rseqc, sam)
from bcbio.broad import BroadRunner, picardrun

import glob
from itertools import product, repeat


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


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


def _group_input_by_cell_type(in_files, delimiter = "_"):
    def _add_entry(d, v):
        base = os.path.basename(v)
        k = base.split(delimiter)[0]
        if "PbN" in k:
            d["PbN"] = d.get("PbN", []) + [v]
        elif "Pb" in k:
            d["Pb"] = d.get("Pb", []) + [v]
        else:
            logger.error("Error grouping by cell type")
            exit(-1)
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
                logger.info("Fixing mate pair information.")
                pairs = combine_pairs(curr_files)
                first = [x[0] for x in pairs]
                second = [x[1] for x in pairs]
                logger.info("Forward: %s" % (first))
                logger.info("Reverse: %s" % (second))
                fixed = view.map(fastq.fix_mate_pairs_with_config,
                                 first, second, [config] * len(first))
                curr_files = list(flatten(fixed))

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
                bamfiles = view.map(sam.sam2bam, tophat_outputs)
                bamsort = view.map(sam.bamsort, bamfiles)
                view.map(sam.bamindex, bamsort)
                final_bamfiles = bamsort
                curr_files = tophat_outputs

            if stage == "htseq-count":
                _emit_stage_message(stage, curr_files)
                htseq_config = _get_stage_config(config, stage)
                htseq_args = zip(*product(curr_files, [config], [stage]))
                htseq_outputs = view.map(htseq_count.run_with_config,
                                         *htseq_args)
                htseq_outdict[condition] = htseq_outputs

            if stage == "coverage":
                logger.info("Calculating RNASeq metrics on %s." % (curr_files))
                nrun = len(curr_files)
                ref = prepare_ref_file(config["stage"][stage]["ref"], config)
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

            if stage == "rseqc":
                _emit_stage_message(stage, curr_files)
                rseqc_config = _get_stage_config(config, stage)
                rseq_args = zip(*product(curr_files, [config]))
                view.map(rseqc.bam_stat, *rseq_args)
                view.map(rseqc.genebody_coverage, *rseq_args)
                view.map(rseqc.junction_annotation, *rseq_args)
                view.map(rseqc.junction_saturation, *rseq_args)
                RPKM_args = zip(*product(final_bamfiles, [config]))
                RPKM_count_out = view.map(rseqc.RPKM_count, *RPKM_args)
                RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
                                            RPKM_count_out)
                """
                                annotate_args = zip(*product(RPKM_count_fixed,
                                             ["gene_id"],
                                             ["ensembl_gene_id"],
                                             ["human"]))
                view.map(annotate.annotate_table_with_biomart,
                         *annotate_args)
                         """
                view.map(rseqc.RPKM_saturation, *rseq_args)
                curr_files = tophat_outputs

    # combine htseq-count files and run deseq on them
    conditions, htseq_files = dict_to_vectors(htseq_outdict)
    deseq_config = _get_stage_config(config, "deseq")
    cell_types = _group_input_by_cell_type(htseq_files)
    for cell_type, files in cell_types.items():
        for comparison in deseq_config["comparisons"]:
            comparison_name = "_vs_".join(comparison)
            deseq_dir = os.path.join(results_dir, "deseq", cell_type,
                                     comparison_name)
            safe_makedir(deseq_dir)
            out_file = os.path.join(deseq_dir, comparison_name + ".counts.txt")
            files_by_condition = _group_input_by_condition(files)
            _emit_stage_message("deseq", files_by_condition)
            c, f = dict_to_vectors(files_by_condition)
            combined_out = htseq_count.combine_counts(f,
                                                      None,
                                                      out_file)
            deseq_out = os.path.join(deseq_dir, comparison_name)
            logger.info("Running deseq on %s with conditions %s "
                        "and writing ot %s" % (combined_out,
                                               conditions,
                                               deseq_out))
            deseq_out = view.map(deseq.run, [combined_out], [c], [deseq_out])
            annotate.annotate_table_with_biomart(deseq_out[0],
                                                 "id",
                                                 "ensembl_gene_id",
                                                 "human")
            #annotated_file = view.map(annotate.annotate_table_with_biomart,
            #                          [deseq_out],
            #                          ["id"],
            #                          ["ensembl_gene_id"],
            #                          ["human"])

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
