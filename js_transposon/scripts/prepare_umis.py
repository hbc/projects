#!/usr/bin/env python
"""Extract genome reads and UMIs from paired end inputs.

- UMI barcode strategy

 5' end of read 1 - L1

 CGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTA

 5' end of read 2 -- needs to be reverse complemented

 A                 UMI       B
 AGTGGCACAGCAGTTAGGNNNNNNNNNNGGGCTGGTCG

Read looks like

L1 -- genome DNA -- L2B -- UMI -- L2A

Algorithm:
- Run PEAR (http://sco.h-its.org/exelixis/web/software/pear/doc.html)
  to merge paired reads.
- Use read structure to extract genomic DNA and UMI into separate FastQ files.

Usage:
   prepare_umis.py <end1_fastq> <end2_fastq>
"""
import collections
import os
import subprocess
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex

adapters = {"L1": "CGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTA",
            "L2B": "CGACCAGCCC",
            "L2A": "CCTAACTGCTGTGCCACT"}
min_genomic_size = 20

def main(end1_fastq, end2_fastq):
    base = os.path.basename(os.path.commonprefix([end1_fastq, end2_fastq])).replace("_L001_R", "")
    merged_file = _merge_fastq(end1_fastq, end2_fastq, base)
    out_genomic = "%s-genomic.fastq" % base
    out_umi = "%s-umi.fastq" % base
    with open(merged_file) as in_handle:
        with open(out_genomic, "w") as genomic_handle:
            with open(out_umi, "w") as umi_handle:
                for name, seq, qual in FastqGeneralIterator(in_handle):
                    seq_qual = extract_genomic_umi(seq, qual)
                    if seq_qual:
                        genomic_handle.write("@%s\n%s\n+\n%s\n" % (name, seq_qual.genomic_seq, seq_qual.genomic_qual))
                        umi_handle.write("@%s\n%s\n+\n%s\n" % (name, seq_qual.umi_seq, seq_qual.umi_qual))

def extract_genomic_umi(seq, qual):
    SeqQual = collections.namedtuple("SeqQual", "genomic_seq,genomic_qual,umi_seq,umi_qual")
    l1_m = regex.search("(%s){s<=4}" % adapters["L1"], seq)
    l2b_m = regex.search("(%s){s<=1}" % adapters["L2B"], seq)
    l2a_m = regex.search("(%s){s<=2}" % adapters["L2A"], seq)
    if l1_m and l2b_m and l2a_m:
        genomic_seq = seq[l1_m.end():l2b_m.start()]
        umi_seq = seq[l2b_m.end():l2a_m.start()]
        genomic_qual = qual[l1_m.end():l2b_m.start()]
        umi_qual = qual[l2b_m.end():l2a_m.start()]
        if len(genomic_seq) >= min_genomic_size:
            return SeqQual(genomic_seq, genomic_qual, umi_seq, umi_qual)
    return None

def _merge_fastq(fq1, fq2, base):
    """Extract paired end reads and merge using pear.
    """
    out_dir = "merged"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    out_base = os.path.join(out_dir, "%s-merge" % (base))
    out_file = out_base + ".assembled.fastq"
    if not os.path.exists(out_file):
        cmd = "pear -v 20 -n 100 -u 1 -f {fq1} -r {fq2} -o {out_base}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])

