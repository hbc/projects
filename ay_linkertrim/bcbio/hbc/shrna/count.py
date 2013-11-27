"""Prepare count information for raw reads within shRNA target regions.
"""
import os
import csv
import collections

from bcbio.utils import file_exists, safe_makedir

import pybedtools

def overlap_target_counts(bam_file, target_file, config):
    """Overlap BAM alignment file with shRNA targets.
    """
    out_dir = safe_makedir(config["dir"]["counts"])
    out_file = os.path.join(out_dir,
                            "{0}.bed".format(os.path.splitext(os.path.basename(bam_file))[0]))
    if not file_exists(out_file):
        pybedtools.set_tempdir(out_dir)
        bed_read_file = pybedtools.BedTool(bam_file).bam_to_bed()
        counts = pybedtools.BedTool(target_file).intersect(bed_read_file, c=True)
        counts.saveas(out_file)
    return out_file

def combine_counts_by_position(count_files, config):
    """Combine multiple count files into a single tab delimited file.
    """
    exp_names = [x["name"] for x in config["experiments"]]
    assert len(count_files) == len(exp_names)
    out_file = os.path.join(os.path.dirname(count_files[0]),
                            "{0}-counts.tsv".format(config["experiment_name"]))
    if not file_exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(["space", "start", "end", "ensembl"] + exp_names)
            readers = [_count_reader(x) for x in count_files]
            base_reader = readers.pop(0)
            for (space, start, end, ensembl, count1) in base_reader:
                other_counts = [x.next()[-1] for x in readers]
                writer.writerow([space, start, end, ensembl, count1] + other_counts)
    return out_file

def _count_reader(in_file):
    with open(in_file) as in_handle:
        for line in in_handle:
            parts = line.rstrip("\r\n").split("\t")
            yield parts

def filter_counts(count_file, config):
    """Filter targets file by multiple criteria.
    """
    fname = _filter_by_gene_only(count_file, config)
    fname = _filter_by_min_experiments(fname, config)
    fname = _filter_by_gene_targets(fname, config)
    return fname

def _filter_count_with_fn(count_file, out_file, test_fn):
    with open(count_file) as in_handle:
        with open(out_file, "w") as out_handle:
            header = in_handle.next()
            out_handle.write(header)
            for line in in_handle:
                if test_fn(line):
                    out_handle.write(line)

def _filter_by_gene_only(count_file, config):
    """Filter targets to only those that have an annotated gene region.
    """
    out_file = apply("{0}-fgene{1}".format, os.path.splitext(count_file))
    def _has_gene_target(line):
        parts = line.rstrip("\r\n").split("\t")
        gene_id = parts[3]
        return gene_id != "."
    if not file_exists(out_file):
        _filter_count_with_fn(count_file, out_file, _has_gene_target)
    return out_file

def _filter_by_gene_targets(count_file, config):
    """Filter target file to gene regions that have multiple targets.
    """
    out_file = apply("{0}-fmulti{1}".format, os.path.splitext(count_file))
    def _get_gene_ids(line):
        return line.rstrip("\r\n").split("\t")[3].split(";")
    if not file_exists(out_file):
        gene_targets = collections.defaultdict(int)
        with open(count_file) as in_handle:
            in_handle.next() # header
            for line in in_handle:
                for gene_id in _get_gene_ids(line):
                    gene_targets[gene_id] += 1
        def _has_multi_targets(line):
            for gene_id in _get_gene_ids(line):
                if gene_targets[gene_id] > 1:
                    return True
            return False
        _filter_count_with_fn(count_file, out_file, _has_multi_targets)
    return out_file

def _filter_by_min_experiments(count_file, config):
    """Filter a count file by reads in multiple experiments with a minimum count.
    """
    min_experiments = config["algorithm"]["count_min_experiments"]
    min_count = config["algorithm"]["count_min_required"]
    def _has_min_experiments(line):
        parts = line.rstrip("\r\n").split("\t")
        counts = [int(x) for x in parts[4:]]
        pass_counts = [x for x in counts if x >= min_count]
        return len(pass_counts) >= min_experiments
    out_file = apply("{0}-fcount{1}".format, os.path.splitext(count_file))
    if not file_exists(out_file):
        _filter_count_with_fn(count_file, out_file, _has_min_experiments)
    return out_file
