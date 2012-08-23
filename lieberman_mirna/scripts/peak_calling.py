import yaml
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import setup_logging
from rkinf.toolbox import macs
import glob
import os
import sys


def main(config_file):

    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    from rkinf.log import logger
    start_cluster(config)
    from rkinf.cluster import view

    input_dir = config["dir"]["input_dir"]
    results_dir = config["dir"]["results"]
    input_files = glob.glob(os.path.join(input_dir, "*.bam"))

    """ example running with macs
    macs.run_with_config(input_file, config, control_file=None, stage=None)
    """

    curr_files = input_files
    for stage in config["run"]:
        # for now just run macs on all of these files without the control
        # file
        if stage == "macs":
            nfiles = len(curr_files)
            out_files = view.map(macs.run_with_config, curr_files,
                                 [config] * nfiles,
                                 [None] * nfiles,
                                 [stage] * nfiles)
            # just use the peak files going forward
            peak_files = [x[0] for x in out_files]
            curr_files = peak_files

        if stage == "intersect":
            """ 1) loop over the ids in the negative group
            """ for each one pick out the files that match it
            """ combine them into one file
            """ output it as the union
            """ 2) loop over the ids in the positive and test group
            """ find intersections of the ones that match the same id:
            """ intersectBed -wao -bed -f fraction -r -a bed1 -b -bed2
            """ might have to try a range of f

    stop_cluster()


if __name__ == "__main__":
    main(sys.argv[1])
