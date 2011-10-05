"""Distributed processing entry points for linker processing.
"""
import os
import subprocess

from celery.task import task

@task(queue="hbc.trim")
def trim_with_aligner(fname, five_trim, three_trim, config):
    """Use bowtie aligner to trim sequence, splitting into aligned sequences.
    """
    assert config["algorithm"]["aligner"] == "bowtie"
    base = os.path.splitext(os.path.basename(fname))[0]
    out_sam = os.path.join(config["dir"]["trim"], "{base}-{trim}.bowtie".format(
        base=base, trim=five_trim))
    out_unaligned = os.path.join(config["dir"]["fastq"],
                                 "{base}-{trim}no.fastq".format(base=base,
                                                                        trim=five_trim))
    if not os.path.exists(out_sam) and not os.path.exists(out_unaligned):
        cl = [config["program"]["bowtie"], "-v", str(config["algorithm"]["align_mismatches"]),
              "-5", str(five_trim), "-3", str(three_trim),
              "--un", out_unaligned, config["algorithm"]["genome"], fname, out_sam]
        print " ".join(cl)
        subprocess.check_call(cl)
    return [{"aligned": covert_bowtieout_to_fastq(out_sam),
             "unaligned": out_unaligned}]

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
