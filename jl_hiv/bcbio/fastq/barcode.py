"""Extract barcodes from fastq input files.
"""
import os
import subprocess
import glob

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.utils import (tmpfile)

# # Demultiplexing

def demultiplex(input_file, barcodes, tmp_dir, config):
    """Polymorphic high level function based on configuration settings.
    """
    prog = config["program"].get("barcode", "")
    if prog.endswith("sabre"):
        return sabre_demultiplex(input_file, barcodes, tmp_dir, config)
    else:
        return bcbb_demultiplex(input_file, barcodes, tmp_dir, config)

# ## sabre barcode de-multiplexing:
# https://github.com/najoshi/sabre

def _write_sabre_bcfile(barcodes, base_file, out_file):
    base, ext = os.path.splitext(base_file)
    unmatched_file = "%s-%s%s" % (base, "unmatched", ext)
    out_files = []
    with open(out_file, "w") as out_handle:
        for name, barcode in barcodes:
            bc_file = "%s-%s%s" % (base, name, ext)
            out_files.append(bc_file)
            out_handle.write("%s %s\n" % (barcode, bc_file))
    return out_files, unmatched_file

def sabre_demultiplex(input_file, barcodes, tmp_dir, config):
    """Do barcode de-multiplexing using sabre.

    Sabre only appears to trim off the 5' side of the read so
    currently not supported.
    """
    raise NotImplementedError
    with tmpfile(dir=tmp_dir, prefix="sabrebc") as bc_file:
        out_files, unmatched_file = _write_sabre_bcfile(barcodes, input_file, bc_file)
        if not os.path.exists(unmatched_file) and not os.path.exists(out_files[0]):
            cl = [config["program"]["barcode"], "se",
                  "-m", str(config["algorithm"]["barcode_mismatch"]),
                  "-f", input_file,
                  "-b", bc_file,
                  "-u", unmatched_file]
            subprocess.check_call(cl)
    return out_files

# ## Custom de-mulitplexing implemented as part of next gen package
# https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/barcode_sort_trim.py

def _write_bcbb_bcfile(barcodes, out_file):
    with open(out_file, "w") as out_handle:
        for name, barcode in barcodes:
            out_handle.write("%s %s\n" % (name, barcode))

def bcbb_demultiplex(input_file, barcodes, tmp_dir, config):
    ext = ".fastq"
    base_name = os.path.splitext(input_file)[0]
    metrics_file = "%s_bc.metrics" % base_name
    out_base = "%s_--b--_--r--%s" % (base_name, ext)
    if not os.path.exists(metrics_file):
        with tmpfile(dir=tmp_dir, prefix="bc") as bc_file:
            _write_bcbb_bcfile(barcodes, bc_file)
            cl = [config["program"]["barcode"], bc_file,
                  out_base, input_file,
                  "--mismatch=%s" % (config["algorithm"]["barcode_mismatch"]),
                  "--metrics=%s" % metrics_file]
            subprocess.check_call(cl)
    return [f for f in glob.glob("%s_*%s" % (base_name, ext))
            if f.find("unmatched") == -1]

# # Barcode conversion and manipulation

def convert_illumina_oldstyle(in_file):
    """Convert older Illumina barcoding conventions to current usage.
    """
    to_remove = ["s_", "_sequence"]
    out_file = in_file
    for rem in to_remove:
        out_file = out_file.replace(rem, "")
    assert out_file != in_file
    if not os.path.exists(out_file):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for name, seq, qual in FastqGeneralIterator(in_handle):
                    bc = name.split("#")[1].split("/")[0]
                    seq += bc
                    qual += qual[0] * len(bc)
                    out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
    return out_file
