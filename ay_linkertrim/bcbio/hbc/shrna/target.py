"""Identify short RNA target locations based on aligned experimental data.

Group together regions of aligned reads from multiple experiments, creating
target regions used to compare counts between experimental conditions.
"""
import os

from bcbio.utils import file_exists, safe_makedir

import pybedtools

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
