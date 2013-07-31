#!/usr/bin/env python
"""Using comparisons to previous array data, identify samples with potential issues.

Uses problem reports from Illumina plus comparisons to arrays to identify
potential sample mismatches.
"""
import glob
import os
import sets
import sys

import pandas as pd
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    df = pd.read_csv(config["array"]["comparison"])
    print "Missing concordance", len(df[df["concordant"] == 0])
    df_prob = df[df["percent"] < 90.0]
    array_prob = set(str(x) for x in df_prob["sample"])
    idremap = read_remap_file(config["idmapping"])
    samplewell_map = sample_to_illuminawell(config["inputs"], idremap)
    illumina_prob = get_illumina_problems(config["illumina"]["problem_report"],
                                          config["inputs"], idremap)
    df["fail_array"] = df["sample"].map(lambda x: str(x) in array_prob)
    df["fail_illumina"] = df["sample"].map(lambda x: str(x) in illumina_prob)
    df_final = df[df.apply(lambda x: x["fail_array"] or x["fail_illumina"], axis=1)]
    new_cols = list(df_final.columns)
    new_cols.insert(1, "illumina_plate")
    df_final["illumina_plate"] = df["sample"].map(lambda x: samplewell_map.get(str(x), ""))
    df_final = df_final.sort_index(by="illumina_plate")
    df_final = df_final.reindex_axis(new_cols, axis=1)
    out_file = "%s-problem%s" % os.path.splitext(config["array"]["comparison"])
    df_final.to_csv(out_file, index=False)
    print "Problem samples written to ", out_file

def get_illumina_problems(problem_report, inputs, idremap):
    """Patient names for samples reported problematic by Illumina.
    """
    out = set()
    with open(problem_report) as in_handle:
        for line in in_handle:
            plate_well = line.rstrip().split("/")[0]
            sample_name = get_input_sample(plate_well, inputs, idremap)
            if sample_name:
                out.add(sample_name["id"])
            else:
                print "Sample not found in current plates: %s" % plate_well
    return out

def sample_to_illuminawell(fpats, idremap):
    out = {}
    for fpat in fpats:
        for dname in glob.glob(fpat):
            sample = dir_to_sample(dname, idremap)
            out[sample["id"]] = os.path.basename(dname)
    return out

def get_input_sample(illumina_well, fpats, idremap):
    for fpat in fpats:
        for dname in glob.glob(fpat):
            if os.path.isdir(dname) and dname.endswith(illumina_well):
                return dir_to_sample(dname, idremap)

def dir_to_sample(dname, idremap):
    vcf_file = os.path.join(dname, "Variations", "SNPs.vcf")
    with open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                illumina_id = line.split("\t")[-1].replace("_POLY", "").rstrip()
                return {"id": idremap.get(illumina_id), "dir": dname,
                        "illuminaid": illumina_id}
    raise ValueError("Did not find sample information in %s" % vcf_file)

def read_remap_file(in_file):
    out = {}
    with open(in_file) as in_handle:
        in_handle.next() # header
        for line in in_handle:
            patient_id, illumina_id = line.rstrip().split()
            out[illumina_id] = patient_id
    return out

if __name__ == "__main__":
    main(sys.argv[1])
