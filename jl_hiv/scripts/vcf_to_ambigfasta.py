#!/usr/bin/env python
"""Convert a list of VCF files into FASTA representation with ambiguous bases.

Usage:
  vcf_to_ambigfasta.py <out file> [<one> <or more> <vcf files>]
"""
import sys
import os

import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import IUPACData

def main(out_file, vcf_files):
    with open(out_file, "w") as out_handle:
        SeqIO.write((vcf_to_rec(f) for f in vcf_files),
                    out_handle, "fasta")

def convert_to_ambig(vcf_rec, ambig_map):
    """Convert information in VCF record into a potentially ambiguous base.
    """
    bases = [x for x in [vcf_rec.REF] + vcf_rec.ALT if x is not None]
    check_bases = "".join(sorted(bases))
    return ambig_map[check_bases]

def vcf_to_rec(in_file):
    """Convert VCF file into Biopython SeqRecord with ambiguous bases.
    """
    ambig_map = get_ambig_map()
    with open(in_file) as in_handle:
        bases = [convert_to_ambig(x, ambig_map) for x in vcf.Reader(in_handle)]
    name = os.path.splitext(os.path.basename(in_file))[0]
    return SeqRecord(Seq("".join(bases)), name, "", "")

def get_ambig_map():
    out = {}
    ignore = "X"
    for k, v in IUPACData.ambiguous_dna_values.iteritems():
        if k not in ignore:
            v = "".join(sorted(list(v)))
            out[v] = k
    return out

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
