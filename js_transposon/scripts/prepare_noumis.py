#!/usr/bin/env python
"""Extract genome reads from single end no-UMI inputs.

- barcode strategy

 5' end of read 1 - L1

 CGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTA

Read looks like

L1 -- genome DNA -- L2B

L2B is not always present, depending on the length of the genomic insert.

Algorithm:
- Use read structure to trim adapter and extract genomic DNA

Usage:
   prepare_umis.py <end1_fastq>
"""
import collections
import gzip
import os
import subprocess
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex

adapters = {"L1": "CGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTA",
            "L2B": "CGACCAG",
            "L2A": "CCTAACTGCTGTGCCACT"}
min_genomic_size = 20

def main(end1_fastq):
    base = _clean_name(os.path.basename(end1_fastq))
    out_genomic = "%s-genomic.fastq" % base
    with gzip.open(end1_fastq) as in_handle:
        with open(out_genomic, "w") as genomic_handle:
            for name, seq, qual in FastqGeneralIterator(in_handle):
                seq_qual = extract_genomic(seq, qual)
                if seq_qual:
                    genomic_handle.write("@%s\n%s\n+\n%s\n" % (name, seq_qual.genomic_seq, seq_qual.genomic_qual))

def _clean_name(name):
    name = name.replace(".fastq.gz", "")
    want = []
    for part in name.split("_"):
        if not (part in ["R1", "001", "S1"] or part.startswith("L00")):
            want.append(part)
    return "_".join(want)

def extract_genomic(seq, qual):
    SeqQual = collections.namedtuple("SeqQual", "genomic_seq,genomic_qual")
    l1_m = regex.search("(%s){s<=4}" % adapters["L1"], seq)
    l2b_m = regex.search("(%s){s<=1}" % adapters["L2B"], seq)
    if l1_m:
        if l2b_m:
            genomic_seq = seq[l1_m.end():l2b_m.start()]
            genomic_qual = qual[l1_m.end():l2b_m.start()]
        else:
            genomic_seq = seq[l1_m.end():]
            genomic_qual = qual[l1_m.end():]
        if len(genomic_seq) >= min_genomic_size:
            return SeqQual(genomic_seq, genomic_qual)
    return None

if __name__ == "__main__":
    main(*sys.argv[1:])

