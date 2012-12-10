import glob
import os
from bipy.utils import flatten
from bcbio.utils import safe_makedir
import sh

def main():
    dirs = ["results/IMPACT_v2/length_filtered/novoalign",
            "results/IMPACT_v3/novoalign"]
    conditions = ["miR-34", "miR-522", "let-7", "control"]
    glob_string = "*.filt-sort.bam"

    files = list(flatten([glob.glob(os.path.join(x, glob_string)) for x in dirs]))
    out_dir = "results/v2_v3_combined/"
    safe_makedir(out_dir)

    for condition in conditions:
        condition_files = [x for x in files if condition in x]
        out_file = os.path.join(out_dir, condition + "_v2_v3.bam")
        print "Combining %s into %s." % (condition_files, out_file)
        bsub_call = list(flatten(["-q", "hsph", "-o", "out" + condition, "-e", "err" + condition, "samtools", "merge", out_file, condition_files]))
        sh.bsub(bsub_call)


if __name__ == "__main__":
    main()
