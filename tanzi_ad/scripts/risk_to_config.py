"""Add risk calculations to a configuration file.

Usage:
  risk_to_config.py <config file> <fam file>
"""
import sys
import os

import yaml

def main(config_file, fam_file):
    risks = _fam_to_risks(fam_file)
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    out = []
    for sample in config["details"]:
        risk = risks[str(sample["description"])]
        sample["metadata"]["risk"] = risk
        out.append(sample)
    config["details"] = out
    out_file = "%s-risks%s" % os.path.splitext(config_file)
    with open(out_file, "w") as out_handle:
        yaml.safe_dump(config, out_handle, default_flow_style=False, allow_unicode=False)

def _fam_to_risks(in_file):
    risks = {}
    with open(in_file) as in_handle:
        header = in_handle.readline().strip().replace("#", "").split("\t")
        header_is = {}
        for i, x in enumerate(header):
            header_is[x] = i
        for line in in_handle:
            parts = line.strip().split("\t")
            risks[parts[header_is["name"]]] = parts[header_is["risk"]]
    return risks

if __name__ == "__main__":
    main(*sys.argv[1:])
