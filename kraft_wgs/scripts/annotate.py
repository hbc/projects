from rkinf.log import setup_logging
from rkinf.cluster import start_cluster, stop_cluster
import yaml


def main(config_file):
    # load yaml config file
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # setup logging
    setup_logging(config)
    from rkinf.log import logger
    # start cluster
    start_cluster(config)

    # grab input files
    # run snpeffects on each in parallel
    # serially load into gemini database
    stop_cluster()

if __name__ == "__main__":
    main(sys.argv[1])
