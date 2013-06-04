"""
ttGAAAG, 301 Day2 MOI0.1, MLV
ttGACTA, FUW Day2 MOI0.1, HIV
ttCATGA, 151 Lib Day10 neg, MLV
ttCAAAC, FUW Day5 MOI0.1, HIV
ttCAGAT, 151 Lib Day20 pos, MLV
ttAGTAC, \#2 Day2 MOI0.1, HIV

We are expecting it to go like this:

R1: possible_contaminant_sequence_barcode_HIV/MLV_genomic_sequence
R2: possible_contaminant_sequence_barcode_adapter_genomic_sequence

The plan:

read in a record with seqio

- Take first 8 bases from each end
- use Bio.pairwise2 module to align each end to each of the barcodes
Use a large gap penalty, we want no gaps. 1 for correct, 1 for mismatch.
- take the best matching barcode as the correct barcode if it is at
least n-1 of match where n is the barcode length
- if none match then read is unclassified (unclassified)
- if both match the same barcode then the read is that barcode (full_evidence)
- if only one matches a barcode then the read is that barcode (half_evidence)
- if each end matches a different barcode, read is unclassified (ambiguous)

pairwise2.align.localxx("barcode", "sequence")
for a in pairwise2.align.localms("ACCGT", "ACG", 1, -1, 20):
give like -20 for a gap, no gaps
-1 for a mismatch
1 for a correct

"""

from Bio import SeqIO, pairwise2
import sys
import yaml
import unittest
from itertools import izip
import os
from collections import defaultdict

config_file = "test_pairs.yaml"
MATCH_LENGTH = 20
OUT_DIR = "barcoded"
CHUNK_SIZE = 100000

def main(config_file, fq1, fq2):
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    barcodes = config["barcodes"]
    anchors = config["anchors"]

    parsers = izip(SeqIO.parse(fq1, "fastq"), SeqIO.parse(fq2, "fastq"))
    fq1_base, fq1_ext = os.path.splitext(fq1)
    fq2_base, fq2_ext = os.path.splitext(fq2)
    barcode_fq1_reads = defaultdict(list)
    barcode_fq2_reads = defaultdict(list)
    count = 0
    for r1, r2 in parsers:
        count = count + 1
        assert r1.id == r2.id
        barcode = _guess_barcode(r1, r2, barcodes, anchors)
        """
        with open(_out_filename(fq1_base, fq1_ext, barcode), "a") as out_handle:
                SeqIO.write(r1, out_handle, "fastq")
        with open(_out_filename(fq2_base, fq2_ext, barcode), "a") as out_handle:
                SeqIO.write(r2, out_handle, "fastq")
                """
        barcode_fq1_reads[barcode].append(r1)
        barcode_fq2_reads[barcode].append(r2)
        if count > CHUNK_SIZE:
            for bc, reads in barcode_fq1_reads.items():
                with open(_out_filename(fq1_base, fq1_ext, bc), "a") as out_handle:
                    SeqIO.write(reads, out_handle, "fastq")
            for bc, reads in barcode_fq2_reads.items():
                with open(_out_filename(fq2_base, fq2_ext, bc), "a") as out_handle:
                    SeqIO.write(reads, out_handle, "fastq")
            count = 0
            barcode_fq1_reads = defaultdict(list)
            barcode_fq2_reads = defaultdict(list)
    for bc, reads in barcode_fq1_reads.items():
        with open(_out_filename(fq1_base, fq1_ext, bc), "a") as out_handle:
            SeqIO.write(reads, out_handle, "fastq")
    for bc, reads in barcode_fq2_reads.items():
        with open(_out_filename(fq2_base, fq2_ext, bc), "a") as out_handle:
            SeqIO.write(reads, out_handle, "fastq")


def _out_filename(base, ext, barcode):
    if not barcode:
        barcode = "ambiguous"
    out_name = os.path.basename(base) + "_" + barcode + ext
    return os.path.join(OUT_DIR, out_name)

def _guess_barcode(r1, r2, barcodes, anchors):
    r1_matches, r2_matches = _match_barcodes(r1.seq[0:MATCH_LENGTH],
                                            r2.seq[0:MATCH_LENGTH], barcodes, anchors)
    return _consensus_barcode(r1_matches, r2_matches)

def _consensus_barcode(b1, b2):
    """
    returns the consensus barcode if any from two sets of barcodes
    """
    s1 = set(b1)
    s2 = set(b2)
    i = s1.intersection(s2)
    if len(i) == 1:
        return i.pop()
    elif len(s1) == 0 and len(s2) == 1:
        return s2.pop()
    elif len(s1) == 1 and len(s2) == 0:
        return s1.pop()
    else:
        return None

def _anchor_barcode(barcode, anchors):
    bc = barcode[0]
    desc = barcode[1]
    return bc + anchors[desc[1]][0], bc + anchors[desc[1]][1],

def _match_barcodes(r1_seq, r2_seq, barcodes, anchors):
    fw_matches = []
    rv_matches = []
    for barcode in barcodes.items():
        fw_anchor, rv_anchor = _anchor_barcode(barcode, anchors)
        score = pairwise2.align.localms(fw_anchor, str(r1_seq), 1, -1, -20, -10, score_only=1)
        if score >= float(len(fw_anchor)) * 0.90:
                fw_matches.append(barcode[0])
        score = pairwise2.align.localms(rv_anchor, str(r1_seq), 1, -1, -20, -10, score_only=1)
        if score >= float(len(fw_anchor)) * 0.90:
                fw_matches.append(barcode[0])
        score = pairwise2.align.localms(fw_anchor, str(r2_seq), 1, -1, -20, -10, score_only=1)
        if score >= float(len(fw_anchor)) * 0.90:
                rv_matches.append(barcode[0])
        score = pairwise2.align.localms(rv_anchor, str(r2_seq), 1, -1, -20, -10, score_only=1)
        if score >= float(len(rv_anchor)) * 0.90:
                rv_matches.append(barcode[0])
                """
        for a in pairwise2.align.localms(fw_anchor, str(r1_seq), 1, -1, -20, -10, score_only=1):
            if a[2] >= float(len(fw_anchor)) * 0.90:
                fw_matches.append(barcode[0])
        for a in pairwise2.align.localms(rv_anchor, str(r1_seq), 1, -1, -20, -10, score_only=1):
            if a[2] >= float(len(rv_anchor)) * 0.90:
                fw_matches.append(barcode[0])
        for a in pairwise2.align.localms(rv_anchor, str(r2_seq), 1, -1, -20, -10, score_only=1):
            if a[2] >= float(len(rv_anchor)) * 0.90:
                rv_matches.append(barcode[0])
        for a in pairwise2.align.localms(fw_anchor, str(r2_seq), 1, -1, -20, -10, score_only=1):
            if a[2] >= float(len(fw_anchor)) * 0.90:
                rv_matches.append(barcode[0])
                """
    return fw_matches, rv_matches

class TestSplitter(unittest.TestCase):

    CONFIG_FILE = "test_pairs.yaml"
    FQ_1 = "test_first10_1.fq"
    FQ_2 = "test_first10_2.fq"

    def setUp(self):
        with open(self.CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
            self.answers = self.config["answers"]
            self.barcodes = self.config["barcodes"]
            fq1_parser = SeqIO.parse(self.FQ_1, "fastq")
            fq2_parser = SeqIO.parse(self.FQ_2, "fastq")
            self.parsers = izip(fq1_parser, fq2_parser)

    def _check_answer(self, read_id, barcode):
        if read_id not in self.answers:
            return False
        else:
            correct = self.answers[read_id]

        assert barcode == correct[0]

    def test_checker(self):
        for r1, r2 in self.parsers:
            self._check_answer(r1.id, self.answers[r1.id][0])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
