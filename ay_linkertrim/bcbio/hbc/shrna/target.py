"""Identify short RNA target locations based on aligned experimental data.

Group together regions of aligned reads from multiple experiments, creating
target regions used to compare counts between experimental conditions.
"""
import os
import operator
import collections

from bcbio.utils import file_exists, safe_makedir

import pybedtools
from pybedtools.cbedtools import parse_attributes
import rpy2.robjects as rpy2

def annotated_target_file(bam_files, config, out_base="shrna_targets"):
    """Prepare target BED file with annotated features.
    """
    bed_targets = identify_targets(bam_files, config, out_base)
    annotations = retrieve_annotations(config)
    ann_targets = add_annotations(bed_targets, annotations)
    return ann_targets

def identify_targets(bam_files, config, out_base="shrna_targets"):
    """Create BED file of target regions based on input BAM alignments
    """
    work_dir = safe_makedir(config["dir"]["annotation"])
    pybedtools.set_tempdir(work_dir)
    out_file = os.path.join(work_dir, "{0}.bed".format(out_base))
    if not file_exists(out_file):
        pybed_files = [pybedtools.BedTool(x) for x in bam_files]
        bed_files = [x.bam_to_bed() for x in pybed_files]
        combined_bed = reduce(lambda x, y: x.cat(y), bed_files)
        merge_bed = combined_bed.merge(d=config["algorithm"].get("merge_distance", 0))
        merge_bed.saveas(out_file)
    return out_file

def retrieve_annotations(config):
    """Use BioMart to retrieve annotation files to intersect with targets.
    """
    out_dir = safe_makedir(config["dir"]["annotation"])
    dataset = config["algorithm"]["biomart_dataset"]
    data_types = ["transcript", "ncRNA"]
    AnnotationFiles = collections.namedtuple("AnnotationFiles", data_types)
    out_files = apply(AnnotationFiles,
                      [os.path.join(out_dir, "{ds}-{type}.gff3").format(ds=dataset, type=x)
                       for x in data_types])
    if not reduce(operator.and_, [file_exists(x) for x in out_files]):
        rpy2.r.assign("out_dir", out_dir)
        rpy2.r.assign("dataset", dataset)
        rpy2.r('''
        options(warn=-1)
        suppressMessages(library(hbc.shrna))
        tx_file <- downloadMartData("transcript", dataset, out_dir)$file
        nc_file <- downloadMartData("ncRNA", dataset, out_dir)$file 
        ''')
    return out_files

def add_annotations(target_file, annotations):
    """Association annotations with BED file of targets.

    Based on the annotate.py example from pybedtools.
    """
    out_file = apply("{0}-annotated{1}".format, os.path.splitext(target_file))
    pybedtools.set_tempdir(os.path.dirname(out_file))
    if not file_exists(out_file):
        all_ann = pybedtools.BedTool(_merge_gff(annotations))
        with_ann = pybedtools.BedTool(target_file).intersect(all_ann, wao=True)
        with open(with_ann.fn) as in_handle:
            with open(out_file, "w") as out_handle:
                _write_combined_features(in_handle, out_handle)
    return out_file

def _write_combined_features(in_handle, out_handle):
    """Prepare BED file with intersecting feature identifiers combined.
    """
    def _write_combined(pos, ids):
        out_handle.write("{0}\t{1}\n".format("\t".join(pos),
                                             ",".join(ids)))
    last = None
    ids = []
    for parts in (l.split("\t") for l in in_handle):
        cur_id = parts[:3]
        cur_attr = parts[-2]
        if cur_id != last:
            if last is not None:
                _write_combined(last, ids)
            last = cur_id
            ids = []
        if cur_attr != ".":
            attrs = parse_attributes(cur_attr)
            ids.append(attrs["ID"])
    if last:
        _write_combined(last, ids)

def _merge_gff(in_files):
    """Merge GFF files into one large combined file.
    """
    merged_file = "{0}merge.gff3".format(os.path.commonprefix(in_files))
    if not file_exists(merged_file):
        with open(merged_file, "w") as out_handle:
            for i, fname in enumerate(in_files):
                with open(fname) as in_handle:
                    for line in in_handle:
                        if line.startswith("#"):
                            if i == 0:
                                out_handle.write(line)
                        else:
                            out_handle.write(line)
    return merged_file
