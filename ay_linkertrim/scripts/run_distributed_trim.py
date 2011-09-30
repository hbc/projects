#!/usr/bin/env python
"""Run distributed linker trimming.

Usage:
  run_distributed_trim.py <system config> <local_config>
"""
import sys

from bcbio import utils
from bcbio.distributed.manage import run_and_monitor

def main(system_config_file, cur_config_file):
    config = utils.merge_config_files([system_config_file, cur_config_file])
    workers_needed = len(config["files"])
    task_module = "bcbio.hbc.linker.tasks"
    queue = "hbc.trim"
    args = [system_config_file, cur_config_file]
    run_and_monitor(config, system_config_file, args, workers_needed, task_module,
                    queues=queue)

if __name__ == "__main__":
    main(*sys.argv[1:])
