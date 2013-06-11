"""Perform trimming of fastq formatted sequence reads.
"""
from contextlib import nested

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.utils import memoize_outfile

@memoize_outfile("-trim.fastq")
def trim_fastq(in_file, five=0, three=0, out_file=None):
    with nested(open(in_file), open(out_file, "w")) as (in_handle, out_handle):
        for name, seq, qual in FastqGeneralIterator(in_handle):
            trim_seq = seq[five:len(seq)-three]
            trim_qual = qual[five:len(qual)-three]
            out_handle.write("@%s\n%s\n+\n%s\n" % (name, trim_seq, trim_qual))
