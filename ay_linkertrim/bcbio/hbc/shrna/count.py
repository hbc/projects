"""Prepare count information for raw reads within shRNA target regions.
"""
import os
import csv

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
    """Filter a count file by reads in multiple experiments with a minimum count.
    """
    min_experiments = config["algorithm"]["count_min_experiments"]
    min_count = config["algorithm"]["count_min_required"]
    def _check_line(line):
        parts = line.rstrip("\r\n").split("\t")
        counts = [int(x) for x in parts[4:]]
        pass_counts = [x for x in counts if x >= min_count]
        return len(pass_counts) >= min_experiments
    out_file = apply("{0}-filter{1}".format, os.path.splitext(count_file))
    if not file_exists(out_file):
        with open(count_file) as in_handle:
            with open(out_file, "w") as out_handle:
                header = in_handle.next()
                out_handle.write(header)
                for line in in_handle:
                    if _check_line(line):
                        out_handle.write(line)
    return out_file
