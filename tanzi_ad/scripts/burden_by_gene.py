"""
script to output for a GEMINI database:
1) a TFAM file
for each gene in the GEMINI database:
2) a TPED file of all variants mapping to a gene
3) the CADD score of all variants in the gene (scaled and raw)

Some of this code would be useful to be added into GEMINI, especially the
constructed ID code as a formatting option. If you are running this in parallel
on a cluster, make sure that all of the CPU you request are on the same node.
"""
import re
import subprocess
import os
import collections
from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count
from functools import partial

def flatten(l):
    """
    flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/

    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

NULL_GENOTYPES = ["."]

def dump_tped(db, out_dir, gene):
    out_file = os.path.join(out_dir, gene + ".tped")
    if os.path.exists(out_file):
        return out_file
    print "Dumping variants in %s to %s." % (gene, out_file)
    cmd = ("gemini query --format tped -q 'select * from variants where "
           "GENE=\"{gene}\"' {db} > {out_file}").format(**locals())
    subprocess.check_call(cmd, shell=True)
    return out_file

def dump_family(db, out_file):
    if os.path.exists(out_file):
        return out_file
    print "Dumping family to %s." % out_file
    cmd = "gemini dump --tfam {db} > {out_file}".format(**locals())
    subprocess.check_call(cmd, shell=True)
    return out_file

def id_from_line(line):
    split_line = line.strip().split("\t")
    chrom = split_line[0].split("chr")[1]
    start = split_line[1]
    end = split_line[2]
    variant_id = split_line[3]
    ref = split_line[4]
    alt = split_line[5]
    gts = split_line[-1]
    geno = [re.split('\||/', x) for x in gts.split(",")]
    geno = [["0", "0"] if any([y in NULL_GENOTYPES for y in x])
            else x for x in geno]
    genotypes = " ".join(list(flatten(geno)))
   # alleles = "|".join(set(list(flatten(geno))).difference("0"))
    return "{chrom}:{start}-{end}:{ref}|{alt}:{variant_id}".format(**locals())

def reformat_cadd_line(line):
    id_string = id_from_line(line)
    split_line = line.strip().split("\t")
    cadd_raw = split_line[6]
    cadd_scaled = split_line[7]
    new_line = "{id_string}\t{cadd_raw}\t{cadd_scaled}"
    return new_line.format(**locals())

def reformat_exon_line(line):
    id_string = id_from_line(line)
    split_line = line.strip().split("\t")
    is_exon = split_line[6]
    new_line = "{id_string}\t{is_exon}"
    return new_line.format(**locals())

def genewise_cadd_scores(db, out_dir, gene):
    out_file = os.path.join(out_dir, gene + ".cadd")
    if os.path.exists(out_file):
        return out_file
    print "Dumping CADD scores in %s to %s." % (gene, out_file)
    cmd = ("gemini query -q 'select chrom, start, end, variant_id, ref, alt, gts, cadd_raw, cadd_scaled from "
           "variants where gene=\"{gene}\"' {db}").format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(out_file, "w") as out_handle:
        out_handle.write("\t".join(["id", "cadd_raw", "cadd_scaled"]) + "\n")
        for line in p.stdout:
            out_handle.write(reformat_cadd_line(line) + "\n")
    return out_file

def dump_exon(db, out_dir, gene):
    out_file = os.path.join(out_dir, gene + ".exon")
    if os.path.exists(out_file):
        return out_file
    cmd = ("gemini query -q 'select chrom, start, end, variant_id, ref, alt, gts, is_exonic from "
           "variants where gene=\"{gene}\"' {db}").format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(out_file, "w") as out_handle:
        out_handle.write("\t".join(["id", "is_exonic"]) + "\n")
        for line in p.stdout:
            out_handle.write(reformat_exon_line(line) + "\n")
    return out_file

def unique_genes(db, sex_only=False):
    if sex_only:
        cmd = ('gemini query -q "select DISTINCT gene from variants where chrom=\'chrX\' or chrom=\'chrY\'" {db}').format(**locals())
    else:
        cmd = ('gemini query -q "select DISTINCT gene from variants" {db}').format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return(set([x.strip() for x in p.stdout]))

if __name__ == "__main__":
    parser = ArgumentParser("Generate gene-wise variant TPED files and CADD scores.")
    parser.add_argument("--out-dir", default="genewise_burden", help="Output directory")
    parser.add_argument("--sex-only", default=False, action='store_true', 
                        help="Only emit sex chromosomes.")
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use.")
    parser.add_argument("db", help="GEMINI database to query.")
    args = parser.parse_args()

    out_dir = args.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print "Calculating unique genes in the database."
    genes = unique_genes(args.db, args.sex_only)
    print "Found %d unique genes." % len(genes)
    print "Dumping TFAM file to ad.tfam."
    dump_family(args.db, "ad.tfam")

    cores_to_use = min(args.cores, cpu_count() / 2)
    print "Using %d cores." % cores_to_use
    p = Pool(processes=cores_to_use)
    partial_cadd = partial(genewise_cadd_scores, args.db, out_dir)
    partial_exon = partial(dump_exon, args.db, out_dir)
    partial_tped = partial(dump_tped, args.db, out_dir)
    print "Dumping CADD scores."
    p.map(partial_cadd, genes)
    print "Dumping exon status."
    p.map(partial_exon, genes)
    print "Dumping TPED."
    p.map(partial_tped, genes)
    print "Finished."
