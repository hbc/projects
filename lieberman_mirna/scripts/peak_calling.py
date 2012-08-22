import yaml
from rkinf.cluster import start_cluster, stop_cluster
from rkinf.log import setup_logging
from rkinf.toolbox import macs
import glob
import os


def main(config_file):

    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    from rkinf.log import logger
    start_cluster(config)

    from rkinf.cluster import view
    input_dir = config["input_dir"]
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
            out_files = view.map(macs.run, curr_files, [config] * nfiles,
                                 [None] * nfiles,
                                 [stage] * nfiles)

            logger.info(out_files)

    stop_cluster()
