#!/usr/bin/env python
"""Download data from SRA and prepare configuration files for pipeline processing.

Usage:
  download_prep.py <YAML prep configuration file>
"""
import sys, os
import urllib
from contextlib import closing
import collections
import subprocess
import datetime

import yaml

from bcbio.utils import safe_makedir, chdir

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    dl_info = parse_sra_metadata(config["data_url"], config["replicates"])
    details = []
    for i, (exp_name, urls) in enumerate(dl_info.iteritems()):
        final_files = [download_fastq(x, config) for x in urls]
        details.append({"lane": i+1,
                        "description": exp_name,
                        "files": final_files,
                        "analysis": config["analysis"],
                        "genome_build": config["genome_build"]})
    write_output_file(details, config)

def write_output_file(details, config):
    out_dir = config["dir"]["config"]
    safe_makedir(out_dir)
    out_file = os.path.join(out_dir, "run_info.yaml")
    out = {"details": details,
           "fc_name": os.path.basename(config["data_url"]),
           "fc_date": str(datetime.datetime.now().strftime("%y%m%d"))}
    with open(out_file, "w") as out_handle:
        yaml.dump(out, out_handle, allow_unicode=False, default_flow_style=False)

def parse_sra_metadata(data_url, replicates_want):
    """Parse experiment metadata from SRA for experiment names and download URLs
    """
    out = collections.defaultdict(list)
    with closing(urllib.urlopen(data_url)) as in_handle:
        in_handle.next() # header
        for line in in_handle:
            parts = line.rstrip("\r\n").split("\t")
            exp_name, _ = parts[7].split()
            dl_url = parts[-1]
            if exp_name.endswith(tuple(replicates_want)):
                out[exp_name].append(dl_url)
    return dict(out)

def download_fastq(url, config):
    out_dir = os.path.abspath(config["dir"]["fastq"])
    safe_makedir(out_dir)
    final_file = os.path.join(out_dir,
                              os.path.splitext(os.path.basename(url))[0])
    if not os.path.exists(final_file):
        with chdir(out_dir):
            subprocess.check_call(["wget", url])
            subprocess.check_call(["gunzip", os.path.basename(url)])
    return os.path.basename(final_file)

if __name__ == "__main__":
    main(*sys.argv[1:])
