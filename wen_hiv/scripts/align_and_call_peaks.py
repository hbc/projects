import subprocess
import sys
import os
import glob

from bcbio.utils import file_exists, replace_suffix
from cluster_helper.cluster import cluster_view




def align(fq1):
    import os
    from bcbio.utils import file_exists, replace_suffix, append_stem
    import subprocess
    genome = "/n/hsphS10/hsphfs1/chb/biodata/genomes/Hsapiens/GRCh37/bowtie2/GRCh37"
    out_sam = os.path.join("align", os.path.basename(replace_suffix(fq1, ".sam")))
    out_bam = replace_suffix(out_sam, ".bam")
    sorted = append_stem(out_bam, "_sorted")
    sorted_prefix = os.path.splitext(sorted)[0]
    out_index = replace_suffix(sorted, ".bai")
    if not file_exists(out_sam):
        cmd = "bowtie2 -S align/{out_sam} {genome} {fq1}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if not file_exists(out_bam):
        cmd = "samtools view -S {out_sam} -b -o {out_bam}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if not file_exists(sorted):
        cmd = "samtools sort {out_bam} {sorted_prefix}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if not file_exists(out_index):
        cmd = "samtools index {sorted}"
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return sorted

def remove_duplicates(sam_file):
    import os
    import subprocess
    from bcbio.utils import file_exists, replace_suffix, append_stem
    md = "/n/HSPH/local/share/java/picard/MarkDuplicates.jar"
    jvm_opts = "-Xms750m -Xmx2000m"
    out_file = append_stem(sam_file, "_dupemarked")
    stats_file = replace_suffix(append_stem(sam_file, "_stats"), ".txt")
    if file_exists(out_file):
        return out_file
    cmd = ("java {jvm_opts} -jar {md} INPUT={sam_file} "
           "OUTPUT={out_file} METRICS_FILE={stats_file} VALIDATION_STRINGENCY=LENIENT")
    subprocess.check_call(cmd.format(**locals()), shell=True)

    return out_file


def main(data_dir, view):
    fastq_files = list(glob.glob(os.path.join(data_dir, "*_trimmed.fastq")))
    print "Aligning %s." % (fastq_files)
    aligned_files = view.map(align, fastq_files)
    print "Deduplicating %s." % (aligned_files)
    deduped = view.map(remove_duplicates, aligned_files)


if __name__ == "__main__":
    data_dir = sys.argv[1]
    fastq_files = list(glob.glob(os.path.join(data_dir, "*_trimmed.fastq")))
    with cluster_view("lsf", "hsph_bioinfo", len(fastq_files), cores_per_job=1) as view:
        main(data_dir, view)
