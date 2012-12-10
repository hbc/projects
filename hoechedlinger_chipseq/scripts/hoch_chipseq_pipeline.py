"""
scripts to run initial fastq mapping for the hochedlinger chipseq analysis

"""
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
import os
from bipy.utils import (append_stem,  flatten, prepare_ref_file,
                        replace_suffix)
from bipy.toolbox import (fastqc, trim, tophat,
                          htseq_count, deseq, fastq, annotate, rseqc, sam)
from bcbio.broad import BroadRunner, picardrun
import glob
from itertools import product


def input_files_from_dir(in_dir, id_file):

    with open(os.path.join(in_dir, id_file)) as in_handle:
        ids = yaml.load(in_handle)

    sample_names = [x for x in ids]
    samples = [glob.glob(in_dir + "/*_%s.R1.fastq" % (x)) for x in sample_names]
    return list(flatten(samples))


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for thesis pipeline
    in_dir = config["dir"]["data"]
    id_file = config["id_file"]
    curr_files = input_files_from_dir(in_dir, id_file)
    logger.info("Running pipeline on %s." % (curr_files))

    for stage in config["run"]:
        if stage == "fastqc":
            logger.info("Running fastqc on %s." % (curr_files))
            stage_runner = fastqc.FastQCStage(config)
            view.map(stage_runner, curr_files)

        if stage == "cutadapt":
            logger.info("Running cutadapt on %s." % (curr_files))
            stage_runner = trim.Cutadapt(config)
            curr_files = view.map(stage_runner, curr_files)

        if stage == "tophat":
            logger.info("Running tophat on %s." % (curr_files))
            tophat_outputs = view.map(tophat.run_with_config,
                                      curr_files, [None] * len(curr_files),
                                      [config["ref"]] * len(curr_files),
                                      ["tophat"] * len(curr_files),
                                      [config] * len(curr_files))
            bamfiles = view.map(sam.sam2bam, tophat_outputs)
            bamsort = view.map(sam.bamsort, bamfiles)
            view.map(sam.bamindex, bamsort)
            final_bamfiles = bamsort
            curr_files = tophat_outputs

        if stage == "coverage":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            nrun = len(curr_files)
            ref = prepare_ref_file(config["stage"][stage]["ref"], config)
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

        if stage == "rseqc":
            logger.info("Running rseqc on %s." % (curr_files))
            rseq_args = zip(*product(curr_files, [config]))
            view.map(rseqc.bam_stat, *rseq_args)
            view.map(rseqc.genebody_coverage, *rseq_args)
            view.map(rseqc.junction_annotation, *rseq_args)
            view.map(rseqc.junction_saturation, *rseq_args)
            RPKM_args = zip(*product(final_bamfiles, [config]))
            RPKM_count_out = view.map(rseqc.RPKM_count, *RPKM_args)
            RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
                                        RPKM_count_out)
            annotate_args = zip(*product(RPKM_count_fixed,
                                         ["gene_id"],
                                         ["ensembl_gene_id"],
                                         ["human"]))
            view.map(annotate.annotate_table_with_biomart,
                     *annotate_args)
            view.map(rseqc.RPKM_saturation, *rseq_args)
            curr_files = tophat_outputs

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
