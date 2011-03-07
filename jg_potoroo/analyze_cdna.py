#!/usr/bin/env python
"""Examine ABI cDNA reads for homology to Marsupials.
"""
import os
import sys
import csv
import glob
import codecs
import itertools
import subprocess

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import rpy2.robjects as rpy2

from bcbio import utils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    _setup_dirs(config)
    clean_cdna = clean_vector(config["cdna_reads"], config)
    cluster_cdna = cluster_seqs(clean_cdna, config)
    for ref in config["ref"]:
        cmp_file = process_ref(cluster_cdna, ref, config)
        annotated_refs(cmp_file, ref, config)

# === Add statistics and descriptions to output file ===

@utils.memoize_outfile("-full.csv")
def annotated_refs(in_file, ref, config, out_file):
    """Use BioMart at Ensembl to add descriptions to each row.
    """
    rpy2.r.assign('in.file', in_file)
    rpy2.r.assign('out.file', out_file)
    rpy2.r.assign("org", ref["ensembl_name"])
    rpy2.r('''
    library(biomaRt)
    options(stringsAsFactors=FALSE)
    in.tbl <- read.csv(in.file, header=TRUE)
    in.tbl$pctsimilar <- in.tbl$hitidentities / in.tbl$hitlength
    print(summary(in.tbl))
    print(sum(in.tbl$hit == ""))

    txs <- unique(in.tbl$hit)
    mart <- useMart("ensembl", dataset=org)
    attrs <- c("ensembl_transcript_id", "description")
    filters <- c("ensembl_transcript_id")
    mart.result <- getBM(attributes=attrs, filters=filters, values=txs, mart=mart)
    names(mart.result) <- c("hit", "description")
    final <- merge(in.tbl, mart.result, by="hit", all.x=TRUE)
    final.sort <- final[order(final$query),]
    print(head(final.sort))
    write.csv(final.sort, out.file, row.names=FALSE, na="")
    ''')

# === Cluster and assemble sequences using wcd and cap3 ===

@utils.memoize_outfile("-cluster.fa")
def cluster_seqs(in_file, config, out_file):
    wcd_out = run_wcd(in_file)
    with open(out_file, "w") as out_handle:
        recs = (r for r in assemble_clusters(in_file, wcd_out, config)
                if len(r.seq) > config["algorithm"]["min_vec_bases"])
        SeqIO.write(recs, out_handle, "fasta")

def assemble_clusters(in_file, wcd_out, config):
    """Provide assembled FASTA records based on wcd clustering.
    """
    rec_find = FastaNumToRec(in_file)
    with open(wcd_out) as wcd_handle:
        for line in wcd_handle:
            nums = [int(n) for n in line.rstrip()[:-1].split()]
            if len(nums) == 1:
                yield rec_find[nums[0]]
            else:
                with utils.tmpfile(prefix="incap3", dir=config["dir"]["work"]) as input_file:
                    with open(input_file, "w") as input_handle:
                        SeqIO.write((rec_find.shortname_rec(n) for n in nums),
                                    input_handle, "fasta")
                    yield cap3_assemble(input_file, config)
                    for fname in glob.glob("%s.cap.*" % input_file):
                        os.remove(fname)

def cap3_assemble(in_file, config):
    """Assemble a FASTA file of clustered sequences with CAP3.
    """
    with utils.tmpfile(prefix="outcap3", dir=config["dir"]["work"]) as cap3_file:
        with open(cap3_file, "w") as out_handle:
            cl = ["cap3", in_file]
            subprocess.check_call(cl, stdout=out_handle)
        seqs = []
        with open(cap3_file) as in_handle:
            for line in in_handle:
                if line.startswith("consensus"):
                    seqs.append(line.rstrip().split()[-1])
    with open(in_file) as in_handle:
        names = []
        for rec in SeqIO.parse(in_handle, "fasta"):
            names.append(rec.id)
    return _make_seqrec("-".join(names), "".join(seqs))

def _make_seqrec(name, seq):
    return SeqRecord(Seq(seq), name, "", "")

@utils.memoize_outfile("-cluster.wcd")
def run_wcd(in_file, out_file):
    cl = ["wcd", "-N", str(config["algorithm"]["cores"]), "-o", out_file,
          "-c", in_file]
    subprocess.check_call(cl)

class FastaNumToRec:
    """Retrieve FASTA records based on index in a file: wcd id -> record.
    """
    def __init__(self, in_file):
        self._name_to_rec = SeqIO.index(in_file, "fasta")
        self._num_to_name = {}
        with open(in_file) as in_handle:
            for i, rec in enumerate(SeqIO.parse(in_handle, "fasta")):
                self._num_to_name[i] = rec.id
    def __getitem__(self, num):
        return self._name_to_rec[self._num_to_name[num]]
    def shortname_rec(self, num):
        rec = self[num]
        rec.id = _parse_query_name(rec.description, True)
        rec.description = ""
        return rec

# === Remove vector ===

@utils.memoize_outfile("-novector.fa")
def clean_vector(in_file, config, out_file):
    """Retrieve FASTA sequence without vector contamination.
    """
    vec_info = process_ref(in_file, config["vector"], config)
    to_clean = []
    with open(vec_info) as in_handle:
        reader = csv.reader(in_handle)
        reader.next()
        for (query, _, _, length, _) in reader:
            if length and int(length) > config["algorithm"]["min_vec_bases"]:
               to_clean.append(query)
    to_clean = tuple(to_clean)
    with open(in_file) as in_handle:
        recs = (r for r in SeqIO.parse(in_handle, "fasta")
                if not r.description.startswith(to_clean))
        with open(out_file, "w") as out_handle:
            SeqIO.write(recs, out_handle, "fasta")

## === Align FASTA to reference with BLAST ===

def process_ref(in_file, ref, config):
    ref_file = ref.get("file", ref.get("cdna", None))
    blast_db = prepare_blast_db(ref_file, "nucl")
    out_file = "%s-%s.tsv" % (os.path.splitext(in_file)[0], ref["name"])
    if not os.path.exists(out_file):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["query", "length", "hit", "hitlength", "hitidentities"])
                with utils.cpmap(config["algorithm"]["cores"]) as cpmap:
                    results = cpmap(blast_homology,
                                    ((rec, blast_db, config["dir"]["work"])
                                     for rec in SeqIO.parse(in_handle, "fasta")))
                    for info in results:
                        writer.writerow(info)
    return out_file

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

def _parse_query_name(title, short=False):
    """From PtK2_plate-001_A04.trimmed.seq (Quality-trimmed) to PtK2_plate-001_A04
    """
    want = title.split()[0]
    name = want.split(".")[0]
    if short:
        name = name.split("-")[-1]
    return name

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
