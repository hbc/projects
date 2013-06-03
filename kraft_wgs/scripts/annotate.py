"""
use brad's script to fix the illumina VCF file headers and combine them.
then annotated the combined vcf file with snpeff.
then load each of these into a gemini database

"""
from bipy.log import setup_logging
from bipy.cluster import start_cluster, stop_cluster
from bcbio.variation import effects, genotype
from bcbio.distributed.transaction import file_transaction
import yaml
import sys
from bipy.pipeline.stages import STAGE_LOOKUP
from bcbio.utils import file_exists
from bipy.utils import flatten, append_stem
import os
import csv
import sh
from itertools import groupby


def main(config_file):
    # load yaml config file
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # setup logging
    setup_logging(config)
    from bipy.log import logger
    # start cluster
    start_cluster(config)
    from bipy.cluster import view

    found = sh.find(config["dir"]["data"], "-name", "Variations")
    var_dirs = [str(x).strip() for x in found]
    logger.info("Var_dirs: %s" % (var_dirs))
    in_dirs = map(os.path.dirname, var_dirs)
    logger.info("in_dirs: %s" % (in_dirs))
    # XXX for testing only load 3
    #curr_files = in_dirs[0:5]
    curr_files = in_dirs


    # run the illumina fixer
    logger.info("Running illumina fixer on %s." % (curr_files))
    illf_class = STAGE_LOOKUP.get("illumina_fixer")
    illf = illf_class(config)
    curr_files = view.map(illf, curr_files)

    # sort the vcf files
    def sort_vcf(in_file):
        from bipy.utils import append_stem
        from bcbio.distributed.transaction import file_transaction
        from bcbio.utils import file_exists
        import sh

        out_file = append_stem(in_file, "sorted")
        if file_exists(out_file):
            return out_file
        with file_transaction(out_file) as tmp_out_file:
            sh.vcf_sort(in_file, _out=tmp_out_file)
        return out_file


    # combine
    out_file = os.path.join(config["dir"].get("results", "results"),
                            "geminiloader",
                            "all_combined.vcf")
    logger.info("Combining files %s into %s." % (curr_files, out_file))
    if file_exists(out_file):
        curr_files = [out_file]
    else:
        curr_files = [genotype.combine_variant_files(curr_files, out_file,
                                                     config["ref"]["fasta"],
                                                     config)]

    # break the VCF files up by chromosome for speed
    logger.info("Breaking up %s by chromosome." % (curr_files))
    breakvcf_class = STAGE_LOOKUP.get("breakvcf")
    breakvcf = breakvcf_class(config)
    curr_files = view.map(breakvcf, curr_files)

    # run VEP on the separate files in parallel
    logger.info("Running VEP on %s." % (curr_files))
    vep_class = STAGE_LOOKUP.get("vep")
    vep = vep_class(config)
    curr_files = view.map(vep, list(flatten(curr_files)))

    curr_files = filter(file_exists, curr_files)

    # load the files into gemini not in parallel
    # don't run in parallel

    # sort the vcf files
    logger.info("Sorting %s." % (curr_files))
    curr_files = view.map(sort_vcf, curr_files)
    # don't run the rest of this in parallel, so take the cluster down
    stop_cluster()

    out_file = os.path.join(config["dir"].get("results", "results"),
                            "geminiloader",
                            "all_combined.vep.vcf")
    logger.info("Combining files %s into %s." % (curr_files, out_file))
    if file_exists(out_file):
        curr_files = [out_file]
    else:
        curr_files = [genotype.combine_variant_files(curr_files, out_file,
                                                     config["ref"]["fasta"],
                                                     config)]


    logger.info("Loading %s into gemini." % (curr_files))
    gemini_class = STAGE_LOOKUP.get("geminiloader")
    geminiloader = gemini_class(config)
    curr_files = map(geminiloader, curr_files)
    logger.info("Run complete.")

if __name__ == "__main__":
    main(sys.argv[1])
