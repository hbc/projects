import glob
import os

from bcbio import utils

remote_names_1 = set(["tanzi-alz-p2f-g1", "tanzi-alz-p2f-g2", "tanzi-alz-p2f-g3"])
remote_dir_1 = "/home/bradc/data/alz-finished"
remote_dir_2 = "/mnt/lustre/stripe-4M/variantdata"

local_dirs = "/n/hsphS10/hsphfs2/tanzi_recalled"

rsync_cmd = "rsync -vlt bradc@dell-sanger:"

has_cram = set(glob.glob(os.path.join(local_dirs, "*alz-*/final/*/*.cram")))

def is_extra_bam(x):
    return x.endswith(("-sr.bam", "-disc.bam", "-sr.bam.bai", "-disc.bam.bai"))

added = set([])
with open("transfer_cram.sh", "w") as out_handle:
    out_handle.write("#!/bin/bash\n#SBATCH -t 7-00:00:00\n#SBATCH -p regal\n#SBATCH --mem=8000\n")
    out_handle.write("set -e\nset -o pipefail\n")
    for test in (sorted(glob.glob(os.path.join(local_dirs, "*alz-*/final/*/*.cram.crai"))) +
                 sorted(glob.glob(os.path.join(local_dirs, "*alz-*/final/*/*.bam"))) +
                 sorted(glob.glob(os.path.join(local_dirs, "*alz-*/final/*/*.bam.bai")))):
        cur_cram = "%s.cram" % utils.splitext_plus(test.replace(".crai", "").replace(".bai", ""))[0]
        if not cur_cram in has_cram and not is_extra_bam(test):
            if "redo" not in cur_cram and test.find(".bam") < 0:
                alz_name = [x for x in cur_cram.split("/") if "alz-" in x][0]
                remote_dir = remote_dir_1 if alz_name in remote_names_1 else remote_dir_2
                base_name = cur_cram[cur_cram.find(alz_name):]
                remote_full = os.path.join(remote_dir, base_name)
                out_handle.write("%s%s %s\n" % (rsync_cmd, remote_full, cur_cram))
            elif cur_cram not in added:
                print cur_cram, os.path.basename(test)
            added.add(cur_cram)

print len(added)
