""" filters fastq files by minimum and maximum read length and then counts the
number of 3' Gs in a fastq file of reads run this from the project root
directory like this: python scripts/count_3prime_g.py config/qc.yaml """

import stat
import os
import sys
from bipy.toolbox import (tagdust, novoalign, fastqc, bedtools, blastn,
                           sickle, htseq_count)
from bipy.toolbox.fasta import filter_seqio, apply_seqio
from bipy.log import setup_logging, logger
from bipy.cluster import start_cluster, stop_cluster
import pandas as pd
from itertools import takewhile
import yaml
from functools import partial
from bcbio.utils import safe_makedir
from bipy.utils import append_stem, replace_suffix, remove_suffix
from bcbio.broad import BroadRunner, picardrun
from bcbio.ngsalign import tophat
#import shutil
import sh

MAX_READ_LENGTH = 10000

def _get_short_names(input_files):
    return [remove_suffix(os.path.basename(x)) for x in
            input_files]

def write_ratios(bam1, bam2, count_file):
    rbam1 = int(sh.samtools.view("-c", bam1))
    rbam2 = int(sh.samtools.view("-c", bam2))
    header = "\t".join(["total reads",
                        "mre reads",
                        "percentage mre"]) + "\n"
    line = "\t".join([str(rbam1), str(rbam2),
                     str(float(rbam2) / rbam1)]) + "\n"

    with open(count_file, "w") as out_handle:
        out_handle.write(header)
        out_handle.write(line)

def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _combine_fastqc(fastqc_files):
    """ combines a set of metrics from the summary files from
    fastqc and produces plots
    XXX: NOT IMPLEMENTED
    """
    pass


def _remove_self_tagdust(tagdust_config, input_files):
    """ creates a temporary contaminant fasta file for each input file which
    does not have its own sequence in it, for filtering with tagdust
    """
    short_names = map(_short_name, input_files)
    contam = tagdust_config["contaminants"]

    def suffix_not_in_name(seq, suffix):
        return not suffix in seq.id

    # build the predicate to test if the short name of the input file is
    # in the fasta sequence
    predicates = [partial(suffix_not_in_name, suffix=x) for x in short_names]

    # filter the contaminant fasta file, outputting a separate file for
    # each input file
    filtered_files = [filter_seqio(contam, predicate, suffix)
                      for suffix, predicate in zip(short_names, predicates)]

    return filtered_files


def _short_name(name):
    # get a short name from the filenames
    return "-".join(takewhile(lambda x: not "R_" in x,
                              name.split("-")[1:]))

def main(config_file):

    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    start_cluster(config)

    # after the cluster is up, import the view to i
    from bipy.cluster import view
    input_files = config["input"]
    results_dir = config["dir"]["results"]

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    curr_files = input_files

    ## qc steps
    for stage in config["run"]:
        if stage == "fastqc":
            # run the basic fastqc
            logger.info("Running %s on %s" % (stage, str(curr_files)))
            fastqc_config = config["stage"][stage]
            fastqc_outputs = view.map(fastqc.run, curr_files,
                                      [fastqc_config] * len(curr_files),
                                      [config] * len(curr_files))
            # this does nothing for now, not implemented yet
            summary_file = _combine_fastqc(fastqc_outputs)

        if stage == "trim":
            logger.info("Trimming poor quality ends "
                        " from %s" % (str(curr_files)))
            nlen = len(curr_files)
            min_length = str(config["stage"][stage].get("min_length", 20))

            # trim low quality ends of reads
            # do this dirty for now
            out_dir = os.path.join(results_dir, "trimmed")
            safe_makedir(out_dir)
            out_files = [append_stem(os.path.basename(x), "trim") for
                         x in curr_files]
            out_files = [os.path.join(out_dir, x) for x in out_files]
            # XXX remove the magic number of 10 the length of the
            # minimum read to keep
            out_files = view.map(sickle.run, curr_files,
                                 ["se"] * nlen,
                                 ["sanger"] * nlen,
                                 [min_length] * nlen,
                                 out_files)
            curr_files = out_files

        if stage == "tagdust":
            input_files = curr_files
            # remove tags matching the other miRNA tested
            logger.info("Running %s on %s." % (stage, input_files))
            tagdust_config = config["stage"][stage]
            tagdust_outputs = view.map(tagdust.run, input_files,
                                       [tagdust_config] * len(input_files),
                                       [config] * len(input_files))
            curr_files = [x[0] for x in tagdust_outputs]

        if stage == "filter_length":
            # filter out reads below or above a certain length
            filter_config = config["stage"][stage]
            min_length = filter_config.get("min_length", 0)
            max_length = filter_config.get("max_length", MAX_READ_LENGTH)

            # length predicate
            def length_filter(x):
                return min_length < len(x.seq) < max_length

            # filter the input reads based on length
            # parallelizing this doesn't seem to work
            # ipython can't accept closures as an argument to view.map()
            """
            filtered_fastq = view.map(filter_seqio, tagdust_outputs,
                                      [lf] * len(tagdust_outputs),
                                      ["filt"] * len(tagdust_outputs),
                                      ["fastq"] * len(tagdust_outputs))"""
            out_files = [append_stem(os.path.basename(input_file[0]),
                         "filt") for input_file in tagdust_outputs]
            out_dir = os.path.join(config["dir"]["results"],
                                   "length_filtered")
            safe_makedir(out_dir)
            out_files = [os.path.join(out_dir, x) for x in out_files]

            filtered_fastq = [filter_seqio(x[0], length_filter, y, "fastq")
                              for x, y in zip(tagdust_outputs, out_files)]

            curr_files = filtered_fastq

        if stage == "count_ends":
            logger.info("Compiling nucleotide counts at 3' and 5' ends.")
            # count the nucleotide at the end of each read
            def count_ends(x, y):
                """ keeps a running count of an arbitrary set of keys
                during the reduce step """
                x[y] = x.get(y, 0) + 1
                return x

            def get_3prime_end(x):
                return str(x.seq[-1])

            def get_5prime_end(x):
                return str(x.seq[0])

            def output_counts(end_function, count_file):
                # if the count_file already exists, skip
                outdir = os.path.join(config["dir"]["results"], stage)
                safe_makedir(outdir)
                count_file = os.path.join(outdir, count_file)
                if os.path.exists(count_file):
                    return count_file
                # outputs a tab file of the counts at the end
                # of the fastq files kj
                counts = [reduce(count_ends,
                                 apply_seqio(x, end_function, kind="fastq"),
                                 {}) for x in curr_files]
                df = pd.DataFrame(counts,
                                  index=map(_short_name, curr_files))
                df = df.astype(float)
                total = df.sum(axis=1)
                df = df.div(total, axis=0)
                df["total"] = total
                df.to_csv(count_file, sep="\t")

            output_counts(get_3prime_end, "3prime_counts.tsv")
            output_counts(get_5prime_end, "5prime_counts.tsv")

        if stage == "tophat":
            tophat_config = config["stage"][stage]
            logger.info("Running tophat on %s" % (str(curr_files)))
            nlen = len(curr_files)
            pair_file = None
            ref_file = tophat_config["annotation"]
            out_base = os.path.join(results_dir, "mirna")
            align_dir = os.path.join(results_dir, "tophat")
            config = config
            tophat_files = view.map(tophat.align,
                                    curr_files,
                                    [pair_file] * nlen,
                                    [ref_file] * nlen,
                                    [out_base] * nlen,
                                    [align_dir] * nlen,
                                    [config] * nlen)
            curr_files = tophat_files

        if stage == "novoalign":
            logger.info("Running novoalign on %s" % (str(curr_files)))
            # align
            ref = config["genome"]["file"]
            novoalign_config = config["stage"][stage]
            aligned_outputs = view.map(novoalign.run, curr_files,
                                       [ref] * len(curr_files),
                                       [novoalign_config] * len(curr_files),
                                       [config] * len(curr_files))
            # convert sam to bam, sort and index
            picard = BroadRunner(config["program"]["picard"])
            bamfiles = view.map(picardrun.picard_formatconverter,
                                [picard] * len(aligned_outputs),
                                aligned_outputs)
            sorted_bf = view.map(picardrun.picard_sort,
                                 [picard] * len(bamfiles),
                                 bamfiles)
            view.map(picardrun.picard_index, [picard] * len(sorted_bf),
                     sorted_bf)
            # these files are the new starting point for the downstream
            # analyses, so copy them over into the data dir and setting
            # them to read only
            #data_dir = os.path.join(config["dir"]["data"], stage)
            #safe_makedir(data_dir)
            #view.map(shutil.copy, sorted_bf, [data_dir] * len(sorted_bf))
            #new_files = [os.path.join(data_dir, x) for x in
            #             map(os.path.basename, sorted_bf)]
            #[os.chmod(x, stat.S_IREAD | stat.S_IRGRP) for x in new_files]
            # index the bam files for later use
            #view.map(picardrun.picard_index, [picard] * len(new_files),
            #         new_files)

            curr_files = sorted_bf

        if stage == "new_coverage":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            nrun = len(curr_files)
            ref = blastn.prepare_ref_file(config["stage"][stage]["ref"],
                                          config)
            ribo = config["stage"][stage]["ribo"]
            picard = BroadRunner(config["program"]["picard"])
            out_dir = os.path.join(results_dir, "new_coverage")
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
            curr_files = out_files

        if stage == "coverage":
            gtf = blastn.prepare_ref_file(config["annotation"], config)
            logger.info("Calculating coverage of features in %s for %s"
                        % (gtf, str(sorted_bf)))
            out_files = [replace_suffix(x, "counts.bed") for
                         x in sorted_bf]
            out_dir = os.path.join(results_dir, stage)
            safe_makedir(out_dir)
            logger.info(out_files)
            out_files = [os.path.join(out_dir,
                                      os.path.basename(x)) for x in out_files]
            logger.info(out_files)
            view.map(bedtools.count_overlaps, sorted_bf,
                     [gtf] * len(sorted_bf),
                     out_files)

        if stage == "htseq-count":
            nfiles = len(curr_files)
            htseq_config = _get_stage_config(config, stage)
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     aligned_outputs,
                                     [config] * nfiles,
                                     [stage] * nfiles)
            column_names = _get_short_names(input_files)
            logger.info("Column names: %s" % (column_names))
            out_file = os.path.join(config["dir"]["results"], stage,
                                    "combined.counts")
            combined_out = htseq_count.combine_counts(htseq_outputs,
                                                      column_names,
                                                      out_file)
        if stage == "bedtools_intersect":
            bedfiles = config["stage"]["bedtools_intersect"].get("bed", None)
            out_dir = os.path.join(results_dir, stage)
            safe_makedir(out_dir)
            for bedfile in bedfiles:
                bedbase, bedext = os.path.splitext(bedfile)
                out_files = [remove_suffix(x) for x in sorted_bf]
                out_files = [os.path.join(out_dir, os.path.basename(x)) for x in
                             out_files]
                out_files = ["_vs_".join([x, os.path.basename(bedbase)])
                             for x in out_files]
                out_files = [".".join([x, "bam"]) for x in out_files]
                test_out = map(bedtools.intersectbam2bed, sorted_bf,
                               [bedfile] * len(sorted_bf),
                               [False] * len(sorted_bf),
                               out_files)
                count_files = [replace_suffix(x, "stats") for x in
                               out_files]
                map(write_ratios, sorted_bf, out_files, count_files)

        if stage == "piranha":
            piranha_runner = piranha.PiranhaStage(config)
            out_files = view.map(piranha_runner, curr_files)

    stop_cluster()

if __name__ == "__main__":
    main(*sys.argv[1:])
