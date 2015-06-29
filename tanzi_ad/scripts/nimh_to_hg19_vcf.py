#!/usr/bin/env python
"""Convert an hg18 tped and tfam file into an hg19 VCF file for comparison purposes.

Relies on Illumina's hg19 remapping file for positions of hg18 probes in hg19.
"""
import itertools
import os
import shutil
import subprocess
import sys

from bcbio import utils

def main(tped_file, tfam_file, remap_file, ref_2bit):
    o_tfam_file, o_tped_file = _remap_tped_tfam(tped_file, tfam_file, remap_file)
    print "Converting to VCF"
    vcf_file = _convert_to_vcf(o_tfam_file, o_tped_file, ref_2bit)

def _convert_to_vcf(tfam_file, tped_file, ref_2bit):
    cmd = [sys.executable, os.path.join(os.getcwd(), "bin", "plink_to_vcf.py"),
           tfam_file, tped_file, ref_2bit]
    subprocess.check_call(cmd)

def _remap_tped_tfam(tped_file, tfam_file, remap_file):
    out_tped_file = os.path.join(os.getcwd(), "%s-hg19.tped" % os.path.splitext(os.path.basename(tped_file))[0])
    out_tfam_file = "%s.tfam" % os.path.splitext(out_tped_file)[0]
    if not utils.file_exists(out_tped_file):
        print "Reading hg19 positions"
        hg19_positions = _read_remap_file(remap_file)
        print "Writing new tped file"
        with open(tped_file) as in_handle:
            with open(out_tped_file, "w") as out_handle:
                for line in in_handle:
                    chrom, pid, x, position, rest = line.split(" ", 4)
                    new_chrom, new_position = hg19_positions.get(pid.lower(), (None, None))
                    if new_chrom:
                        out_handle.write("%s %s %s %s %s" % (new_chrom, pid, x, new_position, rest))
        shutil.copy(tfam_file, out_tfam_file)
    return out_tped_file, out_tfam_file

def _read_remap_file(in_file):
    """Reads Illumina hg19 position file, returning items on standard chromosomes.
    """
    out = {}
    with open(in_file) as in_handle:
        header = in_handle.readline()
        for line in in_handle:
            name, chrom, position, _ = line.split("\t", 3)
            try:
                chrom_test = int(chrom)
            except ValueError:
                chrom_test = None
            if chrom_test and chrom_test > 0:
                out[name.lower()] = (chrom, position)
    return out


if __name__ == "__main__":
    main(*sys.argv[1:])
