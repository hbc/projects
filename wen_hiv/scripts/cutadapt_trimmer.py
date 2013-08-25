import sys
import yaml
import glob
import os
import difflib
from bcbio.utils import safe_makedir, flatten
import subprocess
from bcbio.distributed.transaction import file_transaction
from bcbio.bam.fastq import (filter_reads_by_length,
                             filter_single_reads_by_length)


def main(config_file, fastq_dir):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    barcode_info = config["barcodes"]
    print "Processing %s." % (fastq_dir)
    in_files = glob.glob(os.path.join(fastq_dir, "*.fastq"))
    print "Found %s in %s. " % (in_files, fastq_dir)
    print "Combining paired-end files, if found."
    pairs = combine_pairs(in_files)
    print "Calulcated pairs: %s." % (pairs)
    out_files = []
    for pair in pairs:
        barcode = _determine_barcode_from_filename(pair[0])
        print "Detected barcode: %s" % barcode
        if barcode not in barcode_info.keys():
            print "barcode %s not found in the YAML file, skipping." % (barcode)
            continue
        print "Sample ID: %s" % (barcode_info[barcode][0])
        type = barcode_info[barcode][1]
        print "Sample type: %s" % (barcode_info[barcode][1])
        to_trim = config["to_trim"][type]
        cutadapt_dir = "cutadapt"
        print ("Trimming off %s and any bases before it from %s."
               % (to_trim[0], pair[0]))
        out_dir = os.path.join(cutadapt_dir, os.path.basename(pair[0]))
        out_files.append(_trim_from_front(pair[0], to_trim[0]))
        if len(pair) > 1:
            print ("Trimming off %s and any bases before it from %s."
                   % (to_trim[1], pair[1]))
            out_files.append(_trim_from_front(pair[1], to_trim[1]))
    out_files = list(flatten(out_files))
    out_files = combine_pairs(out_files)
    for pair in out_files:
        if len(pair) > 1:
            filter_reads_by_length(pair[0], pair[1], "fastq-sanger")
        else:
            filter_single_reads_by_length(pair[0], "fastq-sanger")

def _is_barcode_known(barcode, config):
    barcodes = config["barcodes"].keys()
    return barcode in barcodes


def _determine_barcode_from_filename(fq):
    base, _ = os.path.splitext(fq)
    return base.split("_")[-1]


def _trim_from_front(fq, sequence):
    out_dir = "cutadapt"
    safe_makedir(out_dir)
    base = os.path.basename(os.path.splitext(fq)[0])
    untrimmed = os.path.join(out_dir, base + "_untrimmed.fastq")
    trimmed = os.path.join(out_dir, base + "_trimmed.fastq")
    info = os.path.join(out_dir, base + "_info.txt")
    with file_transaction([untrimmed, trimmed, info]) as tx_files:
        tx_trimmed = tx_files[1]
        tx_untrimmed = tx_files[0]
        tx_info = tx_files[2]
        cmd = ("cutadapt --front={sequence} --overlap=10 --output={tx_trimmed} "
               "--untrimmed={tx_untrimmed} {fq} > {tx_info}")
        print ("Trimming %s with cutadapt using the command "
               "%s." % (fq, cmd.format(**locals())))
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return trimmed

def combine_pairs(input_files):
    """ calls files pairs if they are completely the same except
    for one has _1 and the other has _2 returns a list of tuples
    of pairs or singles """
    PAIR_FILE_IDENTIFIERS = ["1", "2"]

    pairs = []
    used = []
    for in_file in input_files:
        if in_file in used:
            continue
        for comp_file in input_files:
            if comp_file in used:
                continue
            s = difflib.SequenceMatcher(a=in_file, b=comp_file)
            blocks = s.get_matching_blocks()
            # length 3 means on match in the middle of the string
            if len(s.get_matching_blocks()) is not 3:
                continue
            if comp_file[blocks[0][2]] in PAIR_FILE_IDENTIFIERS:
                if comp_file[blocks[0][2] - 1] == "_":
                    used.append(in_file)
                    used.append(comp_file)
                    pairs.append([in_file, comp_file])
                    break
        if in_file not in used:
            pairs.append([in_file])
            used.append(in_file)

    return pairs
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
