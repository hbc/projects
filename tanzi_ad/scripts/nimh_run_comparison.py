#!/usr/bin/env python
"""Run validations between WGS and array files, per sample.
"""
import csv
import os
import pprint
import sys

from pysam import VariantFile

from bcbio import utils
from bcbio.provenance import do

def main(wgs_file, array_file, ref_file):
    wgs_file = os.path.abspath(wgs_file)
    array_file = os.path.abspath(array_file)
    samples = _find_samples([wgs_file, array_file])
    wgs_samples = _find_samples([wgs_file])
    _write_samples_used(wgs_samples, samples)
    out_file = "%s-comparison.csv" % os.path.splitext(wgs_file)[0]
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample", "concordance", "TP", "FP", "FN", "mismatch"])
        for sample in samples:
            work_dir = utils.safe_makedir(os.path.join("work", sample))
            with utils.chdir(work_dir):
                metrics_file = _run_comparison(sample, wgs_file, array_file, ref_file)
                metrics = _parse_metrics(metrics_file)
                total = metrics["tp"] + metrics["fp"] + metrics["fn"] + metrics["mismatch"]
                concordance = "%.2f" % (float(metrics["tp"]) / float(total) * 100.0 if total > 0 else 0.0)
                writer.writerow([sample, concordance, metrics["tp"], metrics["fp"], metrics["fn"],
                                 metrics["mismatch"]])
                print(sample, concordance, metrics)

def _parse_metrics(in_vcf):
    metrics = {"tp": 0, "fp": 0, "fn": 0, "mismatch": 0, "nocall": 0, "extra": 0}
    with VariantFile(in_vcf) as in_bcf:
        for rec in in_bcf:
            array = rec.samples["TRUTH"].allele_indices
            wgs = rec.samples["QUERY"].allele_indices
            if sum(array) >= 0 and sum(wgs) >= 0:
                if array == wgs:
                    metrics["tp"] += 1
                elif sum(array) == 0 and sum(wgs) > 0:
                    metrics["fp"] += 1
                elif sum(wgs) == 0 and sum(array) > 0:
                    metrics["fn"] += 1
                else:
                    metrics["mismatch"] += 1
            elif sum(array) > 0:
                metrics["nocall"] += 1
            elif sum(wgs) > 0:
                metrics["extra"] += 1
    return metrics

def _write_samples_used(wgs_samples, array_samples):
    out_file = "talz-wgs-array-comparison-samples.csv"
    array_samples = set(array_samples)
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample", "inarray"])
        for s in wgs_samples:
            writer.writerow([s, "1" if s in array_samples else "0"])

def _run_comparison(sample, orig_wgs_file, orig_array_file, ref_file):
    wgs_file = _subset_sample(orig_wgs_file, sample)
    array_file = _subset_sample(orig_array_file, sample)

    out_prefix = "%s-cmp" % sample
    out_file = out_prefix + ".vcf.gz"
    if not utils.file_exists(out_file):
        locs = ",".join([str(x) for x in range(1, 23)])
        cmd = [os.path.join(os.path.dirname(sys.executable), "hap.py"),
               "-o", out_prefix, "-r", ref_file, "-V",
               "--no-fixchr-truth", "--no-fixchr-query", "-l", locs,
               ] + [array_file, wgs_file]
        do.run(cmd, "Run hap.py comparison")
    return out_file

def _subset_sample(in_file, sample):
    base, ext = os.path.splitext(os.path.basename(in_file))
    out_file = "%s-%s.vcf.gz" % (base, sample)
    if not utils.file_exists(out_file):
        py_cl = os.path.join(os.path.dirname(sys.executable), "py")
        cmd = ("""bcftools view -O v -s {sample} -i 'TYPE="snp"' {in_file} | """
               "vcfstreamsort | "
               """{py_cl} -x 'x if x.startswith("#") or (x.split("\t")[3].lower() in set("gcat") """
               """ and x.split("\t")[4].lower() in set("gcat")) else None' """
               " | bgzip -c > {out_file}")
        do.run(cmd.format(**locals()), "Subset by sample")
        cmd = "tabix -f -p vcf {out_file}"
        do.run(cmd.format(**locals()), "Tabix index")
    return out_file

def _find_samples(in_files):
    cur = set([])
    for fname in in_files:
        bcf_in = VariantFile(fname)
        samples = set(bcf_in.header.samples)
        if len(cur) == 0:
            cur = samples
        else:
            cur = cur.intersection(samples)
    return sorted(list(cur))

if __name__ == "__main__":
    main(*sys.argv[1:])
