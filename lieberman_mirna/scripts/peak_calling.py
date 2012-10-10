import yaml
from bcbio.utils import file_exists
from bipy.cluster import start_cluster, stop_cluster
from bipy.log import setup_logging
from bipy.toolbox import macs
from bipy.utils import remove_suffix
import glob
import os
import sys
import sh

def _merge_condition(in_files, condition):
    """
    merge all of the bam files from a condition together
    as recomended in the MACS manual
    """
    condition_files = [filename for filename in in_files if
                       condition in filename]
    if not condition_files:
        return None
    condition_filename = os.path.join(os.path.dirname(condition_files[1]),
                                      condition + "_merged.bam")
    sorted_prefix = remove_suffix(condition_filename) + ".sorted"
    sorted_filename = sorted_prefix + ".bam"
    if file_exists(sorted_filename):
        return sorted_filename

    sh.samtools("merge", condition_filename, condition_files)
    sh.samtools("sort", condition_filename, sorted_prefix)
    sh.samtools("index", sorted_filename)
    return sorted_filename


def main(config_file):

    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    from bipy.log import logger
    start_cluster(config)
    from bipy.cluster import view

    input_dir = config["dir"]["input_dir"]
    results_dir = config["dir"]["results"]
    input_files = glob.glob(os.path.join(input_dir, "*.bam"))

    """ example running with macs
    macs.run_with_config(input_file, config, control_file=None, stage=None)
    """

    curr_files = input_files
    # first combine all the negative controls into one file
    negative_control = _merge_condition(input_files, 
                                        config["groups"]["negative"])
    test_files = [_merge_condition(input_files, condition) for
                  condition in config["groups"]["test"]]
    test_files = [x for x in test_files if x]
    curr_files = test_files
    
    for stage in config["run"]:
        # for now just run macs on all of these files without the control
        # file
        if stage == "macs":
            nfiles = len(curr_files)
            out_files = view.map(macs.run_with_config, curr_files,
                                 [config] * nfiles,
                                 [negative_control] * nfiles,
                                 [stage] * nfiles)
            # just use the peak files going forward
            peak_files = [x[0] for x in out_files]
            curr_files = peak_files

    stop_cluster()


if __name__ == "__main__":
    main(sys.argv[1])
