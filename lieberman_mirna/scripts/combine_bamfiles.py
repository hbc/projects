"""
combine bamfiles, compute coverage statistics and run jellyfish to
find all 5mers-9mers

"""

import glob
import os
from bipy.utils import flatten, append_stem, remove_suffix, replace_suffix
from bipy.toolbox import blastn
from bcbio.utils import safe_makedir, file_exists
from bcbio.broad import BroadRunner, picardrun
import sh
import yaml
import sys


def main(config_file):
    if config_file:
        with open(config_file) as in_handle:
            config = yaml.load(in_handle)

    dirs = config["in_dir"]
    conditions = config["conditions"]
    glob_string = config["glob_string"]

    files = list(flatten([glob.glob(os.path.join(x, glob_string)) for x in dirs]))
    out_dir = config["dir"]["results"]
    safe_makedir(out_dir)

    curr_files = []
    for condition in conditions:
        condition_files = [x for x in files if condition in x]
        out_file = os.path.join(out_dir, condition + "_v2_v3.bam")
        print "Combining %s into %s." % (condition_files, out_file)
        sh.samtools.merge(list(flatten([out_file, condition_files])))
        #        bsub_call = list(flatten(["-q", "hsph", "-o", "out" + condition, "-e", "err" + condition, "samtools", "merge", out_file, condition_files]))
        #sh.bsub(bsub_call)
        sorted_prefix = remove_suffix(out_file) + ".sorted"
        sorted_file = sorted_prefix + ".bam"
        sh.samtools.sort(out_file, sorted_prefix)
        sh.samtools.index(sorted_file)
        mapped_file = append_stem(sorted_file, "mapped")
        sh.samtools.view(sorted_file, F=4, b=True, o=mapped_file)
        sh.samtools.index(mapped_file)

        # find the reads that don't intersect with the rrna
        in_file = mapped_file
        out_file = os.path.join(out_dir, condition + "_noribo" + "_v2_v3.bam")
        ribo = config["ribo"]
        print "Filtering %s for rRNA in %s into %s." % (in_file, ribo, out_file)
        sh.bedtools.intersect("-abam", in_file, "-v", "-b", ribo, _out=out_file)
        filtered_file = out_file

        print "Calculating RNASeq metrics on %s." % (out_file)
        in_file = out_file
        ref = blastn.prepare_ref_file(config["stage"]["new_coverage"]["ref"],
                                      config)
        ribo = config["stage"]["new_coverage"]["ribo"]
        picard = BroadRunner(config["program"]["picard"])
        out_dir = os.path.join(config["dir"]["results"], "new_coverage")
        safe_makedir(out_dir)
        out_file = replace_suffix(os.path.basename(in_file), "metrics")
        out_file = os.path.join(out_dir, out_file)
        metrics_file = picardrun.picard_rnaseq_metrics(picard, in_file, ref,
                                                       ribo, out_file)

        jelly_dir = os.path.join(config["dir"]["results"], "jellyfish")
        safe_makedir(jelly_dir)
        # convert the filtered file to fastq for jellyfish counting
        fastq_file = os.path.join(jelly_dir,
                                  os.path.basename(replace_suffix(filtered_file,
                                                                  "fastq")))
        sh.bam2fastx(filtered_file, fastq=True, _out=fastq_file)
        for mer in config["stage"]["jellyfish"]["mer_lengths"]:
            base, _ = os.path.splitext(os.path.basename(fastq_file))
            out_prefix = base + "_%dmer" % (mer)
            out_file = os.path.join(jelly_dir, out_prefix)
            if not file_exists(out_file):
                sh.jellyfish.count(fastq_file,
                                   config["stage"]["jellyfish"]["options"],
                                   m=mer, o=out_file)

if __name__ == "__main__":
    main(sys.argv[1])
