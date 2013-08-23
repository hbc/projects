import sys
import os
import glob
import difflib

from cluster_helper.cluster import cluster_view


def align(pair):
    import os
    from bcbio.utils import file_exists, replace_suffix, append_stem, safe_makedir
    import subprocess
    safe_makedir("align")
    genome = "/n/hsphS10/hsphfs1/chb/biodata/genomes/Hsapiens/GRCh37/bowtie2/GRCh37"
    out_sam = os.path.join("align", os.path.basename(replace_suffix(pair[0], ".sam")))
    out_bam = replace_suffix(out_sam, ".bam")
    sorted = append_stem(out_bam, "_sorted")
    sorted_prefix = os.path.splitext(sorted)[0]
    out_index = replace_suffix(sorted, ".bai")
    if not file_exists(out_sam):
        if len(pair) == 2:
            fq1, fq2 = pair
            cmd = "bowtie2 -S {out_sam} {genome} -1 {fq1} -2 {fq2}"
        else:
            fq1 = pair[0]
            cmd = "bowtie2 -S {out_sam} {genome} {fq1}"
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


def count_starts(sam_file):
    import pysam
    import os
    from bcbio.utils import file_exists, safe_makedir
    from collections import Counter
    out_dir = os.path.join(os.path.dirname(os.path.dirname(sam_file)), "coverage")
    safe_makedir(out_dir)
    out_file = os.path.basename(os.path.splitext(sam_file)[0] + "_coverage.bed")
    out_file = os.path.join(out_dir, out_file)
    if file_exists(out_file):
        return out_file
    samfile = pysam.Samfile(sam_file, "rb")
    starts = Counter()
    for read in samfile:
        if read.is_read2:
            continue
        strand = "-" if read.is_reverse else "+"
        pos = str(read.pos if strand == "+" else read.aend)
        chrom = samfile.getrname(read.tid)
        starts["\t".join([chrom, pos, pos, strand])] +=1
    samfile.close()
    with open(out_file, "w") as out_handle:
        for loc, count in starts.items():
            out_handle.write("\t".join([loc, str(count)]) + "\n")
    return out_file


def filter_duplicates(sam_file):
    from bcbio.utils import file_exists, append_stem
    out_file = append_stem(sam_file, "_deduped")
    if file_exists(out_file):
        return out_file
    import pysam
    samfile = pysam.Samfile(sam_file, "rb")
    outsam = pysam.Samfile(out_file, "wb", template=samfile)
    for read in samfile:
        if read.is_unmapped or read.is_duplicate:
            continue
        outsam.write(read)
    samfile.close()
    outsam.close()
    return out_file

def mark_duplicates(sam_file):
    import subprocess
    from bcbio.utils import file_exists, replace_suffix, append_stem
    fm = "/n/HSPH/local/share/java/picard/FixMateInformation.jar"
    md = "/n/HSPH/local/share/java/picard/MarkDuplicates.jar"
    jvm_opts = "-Xms750m -Xmx2000m"
    mate_fixed_file = append_stem(sam_file, "_matefixed")
    if not file_exists(mate_fixed_file):
        cmd = ("java {jvm_opts} -jar {fm} INPUT={sam_file} "
               "OUTPUT={mate_fixed_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    sam_file = mate_fixed_file
    out_file = append_stem(sam_file, "_dupemarked")
    stats_file = replace_suffix(append_stem(sam_file, "_stats"), ".txt")
    if not file_exists(out_file):
        cmd = ("java {jvm_opts} -jar {md} INPUT={sam_file} "
               "OUTPUT={out_file} METRICS_FILE={stats_file} "
               "VALIDATION_STRINGENCY=LENIENT")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file


def compile_duplicate_stats(sam_file):
    import subprocess
    import os
    from bcbio.utils import file_exists
    cmd = "samtools view -F 1024 {sam_file} | wc -l"
    duplicates = subprocess.check_call(cmd.format(**locals()), shell=True)
    cmd = "samtools view {sam_file} | wc -l"
    total = subprocess.check_call(cmd.format(**locals()), shell=True)
    out_file = os.path.splitext(sam_file)[0] + ".stats.txt"
    if file_exists(out_file):
        return out_file
    with open(out_file, "w") as out_handle:
        out_handle.write("Total reads: %s.\n" % (total))
        out_handle.write("Duplicates: %s.\n" % (duplicates))
        out_handle.write("Percentage duplicates: %f.\n"
                         % (float(total) / float(duplicateS)))
    return out_file



def compute_coverage(sam_file):
    import subprocess
    import os
    from bcbio.utils import file_exists, safe_makedir
    out_dir = os.path.join(os.path.dirname(os.path.dirname(sam_file)), "coverage")
    safe_makedir(out_dir)
    out_file = os.path.basename(os.path.splitext(sam_file)[0] + "_coverage.bed")
    out_file = os.path.join(out_dir, out_file)
    if file_exists(out_file):
        return out_file
    cmd = "bedtools genomecov -ibam {sam_file} -bg > {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

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


def main(data_dir, view):
    fastq_files = list(glob.glob(os.path.join(data_dir, "*_trimmed.fixed.fastq")))
    fastq_files = combine_pairs(fastq_files)
    print "Aligning %s." % (fastq_files)
    aligned_files = view.map(align, fastq_files)
    print "Marking duplicates in %s." % (aligned_files)
    marked = view.map(mark_duplicates, aligned_files)
    print "Filtering duplicates and unmapped reads in %s." % (marked)
    deduped = view.map(filter_duplicates, marked)
    #compute_coverage(deduped)
    print "Computing start sites of %s." % (deduped)
    starts = view.map(count_starts, deduped)


if __name__ == "__main__":
    data_dir = sys.argv[1]
    print combine_pairs(list(glob.glob(os.path.join(data_dir, "*_trimmed.fixed.fastq"))))
    fastq_files = combine_pairs(list(glob.glob(os.path.join(data_dir,
                                                            "*_trimmed.fixed.fastq"))))
    with cluster_view("lsf", "hsph_bioinfo", len(fastq_files), cores_per_job=1) as view:
        main(data_dir, view)
