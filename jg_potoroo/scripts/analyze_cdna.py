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
from contextlib import nested

import yaml
import numpy
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
    utils.create_dirs(config)
    clean_cdna = clean_vector(config["cdna_reads"], config)
    cluster_cdna = cluster_seqs(clean_cdna, config)
    tl_stats = translation_stats(cluster_cdna, config)
    for ref in config["ref"]:
        cmp_file = process_ref(cluster_cdna, ref, config)
        display_tl_results(ref, cmp_file, tl_stats, config)
        write_tl_nohit(cluster_cdna, cmp_file, tl_stats, ref, config)
        annotated_refs(cmp_file, ref, config)

# === Add statistics and descriptions to output file ===

@utils.memoize_outfile("-full.csv")
def annotated_refs(in_file, ref, config, out_file):
    """Use BioMart at Ensembl to add descriptions to each row.
    """
    rpy2.r.assign('in.file', in_file)
    rpy2.r.assign('out.file', out_file)
    rpy2.r.assign("org", ref.get("ensembl_name", ""))
    rpy2.r('''
    library(biomaRt)
    options(stringsAsFactors=FALSE)
    in.tbl <- read.csv(in.file, header=TRUE)
    in.tbl$pctsimilar <- in.tbl$hitidentities / in.tbl$hitlength
    print(summary(in.tbl))
    nohits <- sum(in.tbl$hit == "")
    total <- length(in.tbl$hit)
    print(c("No hits", nohits, "Percent hit", (total - nohits) / total))

    if (org != "") {
        txs <- unique(in.tbl$hit)
        mart <- useMart("ensembl", dataset=org)
        attrs <- c("ensembl_transcript_id", "embl", "description")
        filters <- c("ensembl_transcript_id")
        mart.result <- getBM(attributes=attrs, filters=filters, values=txs, mart=mart)
        names(mart.result) <- c("hit", "genbank.id", "description")
        final <- merge(in.tbl, mart.result, by="hit", all.x=TRUE)
        final.sort <- final[order(final$query),]
        print(head(final.sort))
        write.csv(final.sort, out.file, row.names=FALSE, na="")
    }
    ''')

# === Check translations for sequences ===

def write_tl_nohit(in_file, cmp_file, tl_file, ref, config):
    out_file = "%s-%s-nohit.fa" % (os.path.splitext(in_file)[0], ref["name"])
    with nested(open(in_file), open(cmp_file), open(tl_file), open(out_file, "w")) as \
             (in_handle, cmp_handle, tl_handle, out_handle):
        cmp_reader = csv.reader(cmp_handle)
        cmp_reader.next()
        tl_reader = csv.reader(tl_handle)
        tl_reader.next()
        for rec in SeqIO.parse(in_handle, "fasta"):
            (_, size, tl_size) = tl_reader.next()
            (_, _, _, hit) = cmp_reader.next()[:4]
            percent_tl = min(float(tl_size) / float(size), 1.0)
            if not hit and percent_tl >= config["algorithm"]["no_hit_tl"]:
                SeqIO.write([rec], out_handle, "fasta")

def display_tl_results(ref, cmp_file, tl_file, config):
    """Display summary of translated regions versus comparisons.
    """
    def above_thresh(vals):
        above = [v for v in vals if v >= config["algorithm"]["no_hit_tl"]]
        return "Above %s %s" % (config["algorithm"]["no_hit_tl"], len(above))
    hit_stats = []
    nohit_stats = []
    with open(cmp_file) as cmp_handle:
        with open(tl_file) as tl_handle:
            cmp_reader = csv.reader(cmp_handle)
            tl_reader = csv.reader(tl_handle)
            cmp_reader.next()
            tl_reader.next()
            for (_, _, _, hit) in (p[:4] for p in cmp_reader):
                (_, size, tl_size) = tl_reader.next()
                percent_tl = min(float(tl_size) / float(size), 1.0)
                if float(size) > config["algorithm"]["min_stat_size"]:
                    if hit:
                        hit_stats.append(percent_tl)
                    else:
                        nohit_stats.append(percent_tl)
    print "Percent in ORF", ref["name"]
    print "Hits", len(hit_stats), above_thresh(hit_stats)
    print rpy2.r.summary(rpy2.FloatVector(hit_stats))
    print "No hits", len(nohit_stats), above_thresh(nohit_stats)
    print rpy2.r.summary(rpy2.FloatVector(nohit_stats))

@utils.memoize_outfile("-tlstats.csv")
def translation_stats(in_file, config, out_file):
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["read", "size", "tlsize"])
            for rec in SeqIO.parse(in_handle, "fasta"):
                tl_length = _longest_frame(rec, config["dir"]["work"])
                writer.writerow([_parse_query_name(rec.id), len(rec.seq), tl_length])
    return out_file

def _longest_frame(rec, work_dir):
    """Find the longest translatable frame using EMBOSS sixpack.
    """
    lengths = []
    with utils.tmpfile(prefix="insix", dir=work_dir) as in_file:
        with utils.tmpfile(prefix="outsix", dir=work_dir) as out_file:
            with open(in_file, "w") as out_handle:
                SeqIO.write([rec], out_handle, "fasta")
            cl = ["sixpack", "-sequence", in_file, "-outseq", out_file,
                  "-outfile", "/dev/null"]
            with open("/dev/null", "w") as out:
                subprocess.check_call(cl, stderr=out)
            with open(out_file) as in_handle:
                for rec in SeqIO.parse(in_handle, "fasta"):
                    lengths.append(len(rec.seq))
    return max(lengths) * 3

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
    return _make_seqrec("-".join(names), "".join(seqs).replace("-", ""))

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

## === Align FASTA to reference sequences ===

def process_ref(in_file, ref, config):
    db_info = {"blast" : (prepare_blast_db, blast_search),
               "blat" : (prepare_blat_db, blat_search)}
    prepare_db, do_search = db_info[ref.get("aligner", "blast")]
    ref_file = prepare_ref_file(ref, config)
    blast_db = prepare_db(ref_file, "nucl")
    out_file = "%s-%s.tsv" % (os.path.splitext(in_file)[0], ref["name"])
    if not os.path.exists(out_file):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["query", "length", "hit", "hitlength", "hitidentities"])
                with utils.cpmap(config["algorithm"]["cores"]) as cpmap:
                    results = cpmap(do_search,
                                    ((rec, blast_db, config["dir"]["work"])
                                     for rec in SeqIO.parse(in_handle, "fasta")))
                    for info in results:
                        writer.writerow(info)
    return out_file

# ==== Alignments with BLAST ====

@utils.map_wrap
def blast_search(rec, blast_db, tmp_dir):
    with utils.tmpfile(prefix="in", dir=tmp_dir) as in_file:
        with open(in_file, "w") as out_handle:
            SeqIO.write([rec], out_handle, "fasta")
        with utils.tmpfile(prefix="out", dir=tmp_dir) as blast_out:
            return _do_blast(in_file, blast_db, blast_out)

def _do_blast(in_file, blast_db, out_file):
    cl = NcbiblastnCommandline(query=in_file, db=blast_db, out=out_file,
                               outfmt=5, num_descriptions=0, num_alignments=1,
                               evalue=0.1)
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

def prepare_blast_db(db_file, dbtype):
    exts = {"prot" : "pin", "nucl": "nin"}
    base, _ = os.path.splitext(db_file)
    if not os.path.exists("%s.%s" % (base, exts[dbtype])):
        cl = ["makeblastdb", "-in", db_file, "-out", base, "-dbtype", dbtype]
        subprocess.check_call(cl)
    return base

# ==== Alignments with BLAT ====

@utils.map_wrap
def blat_search(rec, db, tmp_dir):
    with utils.tmpfile(prefix="inblat", dir=tmp_dir) as in_file:
        with open(in_file, "w") as out_handle:
            SeqIO.write([rec], out_handle, "fasta")
        with utils.tmpfile(prefix="outblat", dir=tmp_dir) as blat_out:
            return _do_blat(in_file, db, blat_out)

def _do_blat(in_file, db, out_file):
    cl = ["blat", db, in_file, out_file]
    subprocess.check_call(cl)
    raise NotImplementedError

def prepare_blat_db(db_file, dbtype):
    assert dbtype == "nucl"
    out_file = "%s.2bit" % (os.path.splitext(db_file)[0])
    if not os.path.exists(out_file):
        cl = ["faToTwoBit", db_file, out_file]
        subprocess.check_call(cl)
    return out_file

# === Retrieve fasta databases based on URLs ===

def prepare_ref_file(ref, config):
    """Get a reference file, either by URL or locally.
    """
    url = ref.get("url", None)
    if url:
        ref_file = _download_ref(url, config["dir"]["ref"])
    else:
        ref_file = ref.get("file", None)
    assert ref_file is not None and os.path.exists(ref_file), ref_file
    return ref_file

def _download_ref(url, ref_dir):
    dl_file = os.path.basename(url)
    ref_file = None
    for supported_ext, extract_cmd in [(".gz", "gunzip")]:
        if dl_file.endswith(supported_ext):
            ref_file = os.path.join(ref_dir, dl_file[:-len(supported_ext)])
            break
    assert ref_file is not None, url
    if not os.path.exists(ref_file):
        with utils.chdir(ref_dir):
            cl = ["wget", url]
            subprocess.check_call(cl)
            cl = [extract_cmd, dl_file]
            subprocess.check_call(cl)
    return ref_file

# === Utility functions ===

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


if __name__ == "__main__":
    main(*sys.argv[1:])
