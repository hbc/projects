#!/usr/bin/env python
"""Examine ABI cDNA reads for homology to Marsupials.
"""
import os
import sys
import csv
import codecs
import itertools
import subprocess

import yaml
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

from bcbio import utils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    _setup_dirs(config)
    for ref in config["ref"]:
        process_ref(ref, config)

def process_ref(ref, config):
    ref_file = ref.get("file", ref.get("cdna", None))
    blast_db = prepare_blast_db(ref_file, "nucl")
    out_file = "%s-%s.tsv" % (os.path.splitext(config["cdna_reads"])[0],
                              ref["name"])
    with open(config["cdna_reads"]) as in_handle:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["query", "length", "hit", "hitlength", "hitidentities"])
            with utils.cpmap(config["algorithm"]["cores"]) as cpmap:
                #recs = [(rec, blast_db, config["dir"]["work"])
                #        for rec in itertools.islice(SeqIO.parse(in_handle, "fasta"), 10)]
                results = cpmap(blast_homology,
                                ((rec, blast_db, config["dir"]["work"])
                                 for rec in SeqIO.parse(in_handle, "fasta")))
                for info in results:
                    writer.writerow(info)

@utils.map_wrap
def blast_homology(rec, blast_db, tmp_dir):
    with utils.tmpfile(prefix="in", dir=tmp_dir) as in_file:
        with open(in_file, "w") as out_handle:
            SeqIO.write([rec], out_handle, "fasta")
        with utils.tmpfile(prefix="out", dir=tmp_dir) as blast_out:
            return _do_blast(in_file, blast_db, blast_out)

def _do_blast(in_file, blast_db, out_file):
    cl = NcbiblastnCommandline(query=in_file, db=blast_db, out=out_file,
                               outfmt=5, num_descriptions=1, num_alignments=1)
    subprocess.check_call(str(cl).split())
    with codecs.open(out_file, encoding="utf-8", errors="replace") as blast_handle:
        rec = NCBIXML.read(blast_handle)
    if len(rec.alignments) > 0:
        hsp = rec.alignments[0].hsps[0]
        align_name = _parse_ensembl_name(rec.alignments[0].title)
        align_length = hsp.align_length
        align_identities = hsp.identities
    else:
        align_name, align_length, align_identities = ("", "", "")
    return (_parse_query_name(rec.query), rec.query_length, align_name,
            align_length, align_identities)

def _parse_query_name(title):
    """From PtK2_plate-001_A04.trimmed.seq (Quality-trimmed) to PtK2_plate-001_A04
    """
    want = title.split()[0]
    return want.split(".")[0]

def _parse_ensembl_name(title):
    for p in title.split():
        if p.startswith("ENS"):
            return p
    return title

def prepare_blast_db(db_file, dbtype):
    exts = {"prot" : "pin", "nucl": "nin"}
    base, _ = os.path.splitext(db_file)
    if not os.path.exists("%s.%s" % (base, exts[dbtype])):
        cl = ["makeblastdb", "-in", db_file, "-out", base, "-dbtype", dbtype]
        subprocess.check_call(cl)
    return base

def _setup_dirs(config):
    for dname in ["work"]:
        d = config["dir"][dname]
        if not os.path.exists(d):
            os.makedirs(d)

if __name__ == "__main__":
    main(*sys.argv[1:])
