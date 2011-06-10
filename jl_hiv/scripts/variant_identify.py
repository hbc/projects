#!/usr/bin/env python
"""Identify variants in barcoded multiplexed HIV samples.

Usage:
  variant_identify.py <YAML config>
"""
import os
import sys
import csv
import glob
import copy
import operator
import subprocess
import collections
import itertools
from contextlib import closing

import yaml
import pysam
import khmer
from Bio.SeqIO.QualityIO import (FastqGeneralIterator,
                                 _phred_to_sanger_quality_str)
from Bio.Seq import Seq

from bcbio.utils import create_dirs, map_wrap, cpmap
from bcbio.broad import BroadRunner
from bcbio.fastq.barcode import demultiplex, convert_illumina_oldstyle
from bcbio.fastq.unique import uniquify_bioplayground
from bcbio.fastq.trim import trim_fastq
from bcbio.fastq.filter import kmer_filter, remove_ns
from bcbio.ngsalign import novoalign
from bcbio.variation import mixed
from bcbio.variation.summarize import print_summary_counts

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    ref_index = novoalign.refindex(config["ref"], kmer_size=13, step_size=1)
    create_dirs(config)
    for cur in config["input"]:
        in_fastq = cur["fastq"]
        if cur.get("old_style_barcodes", False):
            in_fastq = convert_illumina_oldstyle(in_fastq)
        bc_files = demultiplex(in_fastq, [(b["name"], b["seq"]) for b in cur["barcodes"]],
                               config["dir"]["tmp"], config)
        bcs = _assign_bc_files(cur["barcodes"], bc_files)
        with cpmap(config["algorithm"]["cores"]) as cur_map:
            for _ in cur_map(process_fastq, ((bc, ref_index, cur, config, config_file)
                                             for bc in bcs)):
                pass

# ## Processing and analysis pipeline

@map_wrap
def process_fastq(bc, ref_index, cur_config, config, config_file):
    do_realignment = config["algorithm"].get("realignment", "")
    do_kmercorrect = config["algorithm"].get("kmer_correct", "")
    trim_three = config["algorithm"].get("trim_three", "")
    picard = BroadRunner(config["program"]["picard"], config["program"]["gatk"],
                         config["algorithm"]["java_memory"])
    in_file = bc["file"]
    if trim_three:
        in_file = trim_fastq(in_file, three=int(trim_three))
    if do_kmercorrect:
        in_file = remove_ns(in_file)
        in_file = kmer_filter(in_file, do_kmercorrect, config)
    unique_file = uniquify_bioplayground(in_file, config)
    align_sam = novoalign.align(config["dir"]["align"], ref_index, unique_file,
                                qual_format=cur_config.get("format", None))
    align_bam = to_bamsort(align_sam, unique_file, config, config_file)
    if do_realignment == "gatk":
        picard.run_fn("picard_index", align_bam)
        align_bam = picard.run_fn("gatk_realigner", align_bam, config["ref"],
                                  deep_coverage=True)
    picard.run_fn("picard_index", align_bam)
    if config["algorithm"].get("range_params", None):
        call_analyze_multiple(align_bam, bc, in_file, config)
    else:
        call_bases_and_analyze(align_bam, bc, in_file, config)

def call_analyze_multiple(align_bam, bc, in_file, config):
    """Write output from multiple parameter settings in YAML format.

    This sets up an output file with the raw data for post-processing
    analysis.
    """
    call_stats = []
    for cur_params in apply(itertools.product,
                            [config["algorithm"][p]
                             for p in config["algorithm"]["range_params"]]):
        cur_config = copy.deepcopy(config)
        for name, val in zip(config["algorithm"]["range_params"], cur_params):
            cur_config["algorithm"][name] = val
        stats = call_bases_and_analyze(align_bam, bc, in_file, cur_config, memoize=False)
        call_stats.extend(stats)
    out_file = os.path.join(config["dir"]["stats"], "%s.yaml" %
                            os.path.splitext(os.path.basename(in_file))[0])
    with open(out_file, "w") as out_handle:
        yaml.dump(call_stats, out_handle)

def call_bases_and_analyze(align_bam, bc, in_file, config, memoize=True):
    out = []
    if bc.get("call_bases", True):
        call_file, params = position_percent_file(align_bam, in_file, config, memoize)
        if bc.get("control", False):
            out = _print_control_summary(call_file, align_bam, config, params, out)
        else:
            identify_patient_variations(call_file, align_bam, config, params)
    return out

# ## Output stats for assessing reliability

def _print_control_summary(call_file, align_bam, config, params, out):
    for expect in config["expected"]:
        out_info = {"file": align_bam, "region": expect["name"], "calls": []}
        out_info.update(params)
        counts = mixed.compare_files(call_file, expect["file"],
                                     expect["offset"], True)
        _print_expect_info(expect["name"], counts)
        for percent, vals in counts.items():
            vals["percent"] = percent
            out_info["calls"].append(vals)
        out.append(out_info)
        print_summary_counts(out_info)
    return out

def _print_expect_info(name, counts):
    print "** %s" % name
    percents = sorted(counts.keys(), reverse=True)
    print "| Percent | Correct | Wrong (partial) | Wrong |"
    print "|---------+---------+-----------------+-------|"
    for percent in percents:
        print "| % 7s | % 7s | % 15s | % 5s |" % (percent,
                                                  counts[percent].get("correct", 0),
                                                  counts[percent].get("partial", 0),
                                                  counts[percent].get("wrong", 0))

# ## Summary of minor variants in a new patient population

def identify_patient_variations(call_file, align_bam, config, params):
    base_order = ["A", "C", "G", "T"]
    out_file = os.path.join(config["dir"]["calls"],
                            "{0}.csv".format(os.path.splitext(os.path.basename(align_bam))[0]))
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["gene", "pos", "called.base"] + base_order)
        for call_info in _patient_variation_iterator(call_file, config["expected"],
                                                     base_order):
            writer.writerow(call_info)

def _patient_variation_iterator(call_file, expect_configs, base_order):
    for expect in expect_configs:
        print expect["name"]
        for call in mixed.call_expected_iter(call_file, expect["file"],
                                             expect["offset"], True):
            for base, percent in call.called.iteritems():
                if percent is not None and percent < 50:
                    epercent = call.expected[base]
                    if epercent is None or int(epercent) < 50:
                        yield ([expect["name"], call.pos, base] +
                               [call.called[b] for b in base_order])

# ## Variant calling based on k-mer filtering and read analysis

def position_percent_file(align_bam, read_file, config, memoize=True):
    params = {"kmer_size": config["algorithm"]["kmer_size"],
              "kmer": config["algorithm"]["kmer_thresh"],
              "qual": int(config["algorithm"].get("qual_thresh", 0)),
              "align_score": config["algorithm"].get("align_score_thresh", 0),
              "call_thresh": config["algorithm"]["call_thresh"]}
    print align_bam, params
    bases = ["A", "C", "G", "T"]
    out_file = os.path.join(config["dir"]["vrn"], "%s-variations.tsv" %
                            os.path.splitext(os.path.basename(align_bam))[0])
    if not memoize or not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(["space", "pos"] + bases)
            ktable, read_counts = count_kmers_and_reads(read_file, params["kmer_size"])
            for chrom, pos, kmers in positional_kmers(align_bam, params):
                base_percents = {}
                for base, percent in base_kmer_percents(kmers, ktable, read_counts, params):
                    base_percents[base] = "%.1f" % (percent * 100.0)
                writer.writerow([chrom, pos] + [base_percents.get(b, "") for b in bases])
    return out_file, params

def base_kmer_percents(kmers, ktable, read_counts, params):
    """Retrieve percentages of each base call based on k-mer counts.
    """
    param_checker = param_position_checker(params)
    kmer_counts = dict((k, ktable.get(k))
                       for k in list(set(x.kmer for x in kmers)))
    base_counts = collections.defaultdict(int)
    total = float(sum(kmer_counts.values()))
    kmers_w_percents = (k._replace(kmer_percent=kmer_counts[k.kmer] / total if total > 0 else 0)
                        for k in kmers)
    for kmer in filter(param_checker, kmers_w_percents):
        base_counts[kmer.call] += read_counts[kmer.seq]
    pass_total = float(sum(base_counts.values()))
    final = []
    for base, count in base_counts.iteritems():
        percent = count / pass_total
        if percent >= params["call_thresh"]:
            final.append((base, count / pass_total))
    final.sort(key=operator.itemgetter(1), reverse=True)
    return final

def param_position_checker(params):
    """Determine if a position passes all parameter filters.
    """
    def _check_position(kmer):
        if kmer.kmer_percent >= params["kmer"]:
            if kmer.align_score >= params["align_score"]:
                if kmer.qual >= params["qual"]:
                    return True
        return False
    return _check_position

def positional_kmers(in_bam, params):
    """Retrieve informative kmers at each piled up position in an alignment.
    """
    qual_map = {}
    for k, v in _phred_to_sanger_quality_str.iteritems():
        qual_map[v] = k
    with closing(pysam.Samfile(in_bam, 'rb')) as work_bam:
        #for col in (p for p in work_bam.pileup() if p.pos > 862 and p.pos < 891):
        #for col in (p for p in work_bam.pileup() if p.pos > 3280 and p.pos < 3310):
        for col in work_bam.pileup():
            space = work_bam.getrname(col.tid)
            kmers = filter(lambda x: x is not None,
                           [_read_surround_region(r, params, qual_map)
                            for r in col.pileups])
            yield space, col.pos, kmers

def _read_surround_region(read, params, qual_map):
    """Provide context for an aligned read at a particular position.

    Requires a full length kmer at each side so excludes information near
    start and end of reads.
    """
    Kmer = collections.namedtuple('Kmer', ['kmer', 'call', 'seq', 'qual',
                                           'align_score', 'kmer_percent'])
    kmer_size = params["kmer_size"]
    assert kmer_size % 2 == 1, "Need odd kmer size"
    extend = (kmer_size - 1) // 2
    if read.indel == 0:
        seq = read.alignment.seq
        pos = read.qpos
        if (pos >= extend and pos < len(seq) - extend):
            qual = qual_map[read.alignment.qual[read.qpos]]
            align_score = [n for (t, n) in read.alignment.tags if t == "AS"][0]
            call = seq[pos]
            kmer = seq[pos-extend:pos+extend+1]
            # if reverse, return forward original read values for counting
            if read.alignment.is_reverse:
                seq = str(Seq(seq).reverse_complement())
                kmer = str(Seq(kmer).reverse_complement())
            assert len(kmer) == kmer_size, (kmer, seq, pos, len(seq))
            return Kmer(kmer, call, seq, qual, align_score, 0.0)

def count_kmers_and_reads(in_fastq, kmer_size):
    ktable = khmer.new_ktable(kmer_size)
    read_count = collections.defaultdict(int)
    with open(in_fastq) as in_handle:
        i = 0
        for (_, seq, _) in FastqGeneralIterator(in_handle):
            i += 1
            #if i > 1e5: break
            if seq.find("N") == -1:
                ktable.consume(seq)
                read_count[seq] += 1
    return ktable, dict(read_count)

# ## Utility functions

def to_bamsort(sam_file, fastq_file, config, config_file):
    sample_name = os.path.splitext(os.path.basename(sam_file))[0]
    cl = [config["program"]["bamsort"], "--name=%s" % sample_name,
          config_file, sam_file, config["ref"], fastq_file]
    subprocess.check_call(cl)
    return glob.glob("%s*-sort.bam" % os.path.splitext(sam_file)[0])[0]

def _assign_bc_files(bcs, files):
    """Assign file names projected by demultiplexing to the barcode.
    """
    final_bcs = []
    for bc in bcs:
        cur_file = None
        for fname in files:
            if fname.find(bc["name"]) >= 0:
                assert cur_file is None
                cur_file = fname
        assert cur_file is not None
        bc["file"] = cur_file
        final_bcs.append(bc)
    return final_bcs

if __name__ == "__main__":
    main(*sys.argv[1:])
