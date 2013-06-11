#!/usr/bin/env python
"""Summarize transcript counts from Cufflinks in regions of interest.

Usage:
  transcript_analyze.py <prep YAML> <read YAML> <work directory>
"""
import sys
import os
import csv
import glob
import contextlib
import collections
import subprocess

import yaml

from bcbio.utils import safe_makedir
from BCBio import GFF

def main(prep_yaml, read_yaml, work_dir):
    outdir = make_outdir(prep_yaml)
    regions = regions_of_interest(prep_yaml)
    sample_info = []
    tx_order = []
    for sample in tx_files(read_yaml, work_dir):
        tx_counts = collections.defaultdict(list)
        with csv_writer(sample, outdir) as writer:
            for region in regions:
                for tx in parse_tx_in_region(sample.txfile, region):
                    writer(tx)
                    tx_counts[tx["transcript_id"]].append(float(tx["FPKM"]))
                    tx_order.append((min(tx["start"], tx["end"]), tx["transcript_id"]))
        sample_info.append((sample, tx_counts))
    tx_order.sort()
    write_summary(tx_order, sample_info, outdir)

def write_summary(tx_order, sample_info, outdir):
    """Write top level summary of gene expression for each transcript.
    """
    out_file = os.path.join(outdir, "ibd_region_summary.csv")
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["name"] +
                        [x.name.split("_")[-1] for x, _ in sample_info] +
                        ["description"])
        for name, counts, description in summary_details(tx_order, sample_info, outdir):
            writer.writerow([name] + counts + [description])

def summary_details(tx_order, sample_info, outdir):
    """Summarize: genes as rows and samples as columns, with readable names.
    """
    for name, description, tx_ids in ensembl_tx_names(tx_order, outdir):
        sample_counts = []
        for _, counts in sample_info:
            cur_counts = []
            for tx_id in tx_ids:
                cur_counts.extend(counts[tx_id])
            cur_count = max(cur_counts) if len(cur_counts) > 0 else 0.0
            sample_counts.append("%.3f" % cur_count)
        yield name, sample_counts, description

def ensembl_tx_names(tx_order, outdir):
    """Retrieve gene name and descriptions using Ensembl.
    """
    tx_file = os.path.join(outdir, "ensembl_tx_info.csv")
    if not os.path.exists(tx_file):
        run_biomart_script(list(set([t for _, t in tx_order])), tx_file)
    return sort_by_tx_order(tx_file, tx_order)

def sort_by_tx_order(tx_file, tx_order):
    """Sort gene names by locations in transcript order.
    """
    tx_locs = {}
    for loc, tx in tx_order:
        tx_locs[tx] = loc
    out = []
    gene_info = group_tx_by_gene_name(tx_file)
    for name, (tx_ids, desc) in gene_info.iteritems():
        loc = min([tx_locs[x] for x in tx_ids])
        out.append((loc, (name, desc, tx_ids)))
    out.sort()
    return [x for (_, x) in out]

def group_tx_by_gene_name(tx_file):
    """Use Ensembl data to group transcripts into named genes for output.
    """
    by_name = {}
    with open(tx_file) as in_handle:
        reader = csv.reader(in_handle)
        for (_, tx_id, gene_name, gene_description) in reader:
            try:
                cur_ids, _ = by_name[gene_name]
            except KeyError:
                cur_ids = []
            cur_ids.append(tx_id)
            by_name[gene_name] = (cur_ids, gene_description)
    return by_name

def run_biomart_script(tx_ids, out_file):
    """Run R biomaRt script to generate gene names and descriptions.
    """
    cl = ["Rscript", os.path.join(os.path.dirname(__file__), "ensembl_gene_info.R"),
          out_file, ",".join(tx_ids)]
    subprocess.check_call(cl)

@contextlib.contextmanager
def csv_writer(sample, outdir):
    """Prepare writer function to output transcript information as a CSV.
    """
    header = ["gene_id", "FPKM", "transcript_id", "chr", "start", "end"]
    with open(os.path.join(outdir,
                           "{0}.csv".format(sample.name)), "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(header)
        def write_tx(tx):
            writer.writerow([tx[x] for x in header])
        yield write_tx

def parse_tx_in_region(tx_file, region):
    want = ["gene_id", "transcript_id", "FPKM"]
    se_range = set(range(region["start"], region["end"]))
    limit_info = {"gff_id": [region["space"]],
                  "gff_type": ["transcript"]}
    for rec in GFF.parse_simple(tx_file, limit_info=limit_info):
        s, e = rec["location"]
        if s in se_range or e in se_range:
            out = {"chr": rec["rec_id"], "start": s, "end": e}
            for n in want:
                out[n] = rec["quals"][n][0]
            yield out

def tx_files(read_yaml, work_dir):
    """Retrieve references to cufflinks transcript GTF files for each sample.
    """
    Sample = collections.namedtuple("Sample", "name txfile")
    with open(read_yaml) as in_handle:
        config = yaml.load(in_handle)
    out = []
    for c in config["details"]:
        fnames = glob.glob(os.path.join(work_dir,
                                        "{lane}_{date}_{fc}*-cufflinks".format(
                                            lane=c["lane"], date=config["fc_date"],
                                            fc=config["fc_name"]),
                                        "transcripts.gtf"))
        assert len(fnames) == 1
        out.append(Sample(c["description"], fnames[0]))
    return out

def regions_of_interest(prep_yaml):
    """Get chrom, start, end coordinates from config YAML
    """
    with open(prep_yaml) as in_handle:
        config = yaml.load(in_handle)
    return config["regions"]

def make_outdir(prep_yaml):
    with open(prep_yaml) as in_handle:
        config = yaml.load(in_handle)
    return safe_makedir(config["dir"]["tx"])

if __name__ == "__main__":
    main(*sys.argv[1:])
