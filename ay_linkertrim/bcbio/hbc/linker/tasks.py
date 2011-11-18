"""Distributed processing entry points for linker processing.
"""
import os
import subprocess

from Bio import SeqIO

from celery.task import task
from bcbio.distributed.transaction import file_transaction

@task(queue="hbc.trim")
def trim_with_aligner(fname, idx, trim_vals, config):
    """Use bowtie aligner to trim sequence, splitting into aligned sequences.
    """
    # workaround for JSON passing that converts integer keys into strings
    swap = {}
    for k, v in trim_vals.iteritems():
        swap[str(k)] = v
    trim_vals = swap
    five_trim, three_trim = trim_vals[str(_read_size(fname))][idx]
    assert config["algorithm"]["aligner"] == "bowtie"
    base = os.path.splitext(os.path.basename(fname))[0]
    out_sam = os.path.join(config["dir"]["trim"], "{base}-{trim}.bowtie".format(
        base=base, trim=five_trim))
    out_unaligned = os.path.join(config["dir"]["work_fastq"],
                                 "{base}-{trim}no.fastq".format(base=base,
                                                                        trim=five_trim))
    if not os.path.exists(out_sam):
        with file_transaction(out_unaligned, out_sam) as (tx_unaligned, tx_sam):
            cl = [config["program"]["bowtie"], "-v", str(config["algorithm"]["align_mismatches"]),
                  "-5", str(five_trim), "-3", str(three_trim),
                  "--un", tx_unaligned, config["algorithm"]["genome"], fname, tx_sam]
            print " ".join(cl)
            subprocess.check_call(cl)
    return [{"aligned": covert_bowtieout_to_fastq(out_sam),
             "unaligned": out_unaligned if os.path.exists(out_unaligned) else None}]

def _read_size(fname):
    """Read length in the current file, used to define trimming.
    """
    with open(fname) as in_handle:
        rec = SeqIO.parse(fname, "fastq").next()
    return len(rec.seq)

def covert_bowtieout_to_fastq(fname):
    """Convert bowtie aligned output into fastq format.
    """
    out = "{base}.fastq".format(base=os.path.splitext(fname)[0])
    if not (os.path.exists(out) and os.path.getsize(out) > 0):
        with open(fname) as in_handle:
            with open(out, "w") as out_handle:
                for name, _, _, _, seq, qual in (l.split("\t")[:6] for l in in_handle):
                    out_handle.write("@{name}\n{seq}\n+\n{qual}\n".format(
                        name=name, seq=seq, qual=qual))
    return out
