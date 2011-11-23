"""Prepare shRNA targets for differential expression analysis.
"""
import os
import csv
import subprocess

import yaml

from bcbio.utils import file_exists, safe_makedir

def do_comparisons(count_file, config):
    for cmp_info in config["comparisons"]:
        for condition in cmp_info["conditions"]:
            noreplicate_comparison(count_file, condition,
                                   cmp_info["background"], config)

def noreplicate_comparison(count_file, condition, background, config):
    """Prepare a differential expression comparison without replicates.
    """
    cur_config, cur_config_file = _prepare_yaml_config(condition, background, config)
    _prepare_count_file(count_file, cur_config["infile"], condition, background)
    subprocess.check_call(["Rscript", config["program"]["diffexp"], cur_config_file])

def _prepare_count_file(orig_count, new_count, condition, background):
    """Prepare subset count with condition and background of interest.
    """
    with open(orig_count) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        header = reader.next()
        cond_i = header.index(condition)
        back_i = header.index(background)
        with open(new_count, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["shrna.id", condition, background, "accession"])
            for parts in reader:
                target_id = apply("{0}:{1}-{2}".format, parts[:3])
                writer.writerow([target_id, parts[cond_i], parts[back_i], parts[3]])

def _prepare_yaml_config(condition, background, config):
    """Prepare YAML configuration file for input to differential expression.
    """
    tmp_dir = safe_makedir(config["dir"]["tmp"])
    out_dir = safe_makedir(config["dir"]["expression"])
    base = "{exp}-{name}".format(exp=config["experiment_name"], name=condition)
    cur_config = {
        "infile": os.path.join(tmp_dir, "{0}-counts.csv".format(base)),
        "out_base": os.path.join(out_dir, base),
        "id_name": "shrna.id",
        "model": {"condition": [condition, background]}}
    config_file = os.path.join(tmp_dir, "{0}-config.yaml".format(base))
    with open(config_file, "w") as out_handle:
        yaml.dump(cur_config, out_handle)
    return cur_config, config_file
