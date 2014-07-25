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
from multiprocessing import Pool

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

def dump_tped(db, gene, out_dir):
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

def reformat_line(line):
    split_line = line.strip().split("\t")
    chrom = split_line[0]
    start = split_line[1]
    end = split_line[2]
    gts = split_line[5]
    cadd_raw = split_line[3]
    cadd_scaled = split_line[4]
    geno = [re.split('\||/', x) for x in gts.split(",")]
    geno = [["0", "0"] if any([y in NULL_GENOTYPES for y in x])
            else x for x in geno]
    genotypes = " ".join(list(flatten(geno)))
    alleles = "|".join(set(list(flatten(geno))).difference("0"))
    new_line = "{chrom}:{start}-{end}:{alleles}\t{cadd_raw}\t{cadd_scaled}"
    return new_line.format(**locals())

def genewise_cadd_scores(db, gene, out_dir):
    out_file = os.path.join(out_dir, gene + ".cadd")
    if os.path.exists(out_file):
        return out_file
    print "Dumping CADD scores in %s to %s." % (gene, out_file)
    cmd = ("gemini query -q 'select chrom, start, end, cadd_raw, cadd_scaled, gts from "
           "variants where gene=\"{gene}\"' {db}").format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(out_file, "w") as out_handle:
        out_handle.write("\t".join(["id", "cadd_raw", "cadd_scaled"]) + "\n")
        for line in p.stdout:
            out_handle.write(reformat_line(line) + "\n")
    return out_file


def get_cadd_score_runner(db, out_dir):
    def wrapper(gene):
        return genewise_cadd_scores(db, gene, out_dir)
    return wrapper

def get_dump_tped_runner(dn, out_dir):
    def wrapper(gene):
        return dump_tped(db, gene, out_dir)
    return wrapper


def unique_genes(db):
    cmd = ('gemini query -q "select DISTINCT gene from variants" {db}').format(**locals())
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return(set([x.strip() for x in p.stdout]))


if __name__ == "__main__":
    parser = ArgumentParser("Generate gene-wise variant TPED files and CADD scores.")
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use.")
    parser.add_argument("db", help="GEMINI database to query.")
    args = parser.parse_args()

    out_dir = "genewise_burden"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print "Calculating unique genes in the database."
    genes = unique_genes(args.db)
    print "Found %d unique genes." % len(genes)
    print "Dumping TFAM file to ad.tfam."
    dump_family(args.db, "ad.tfam")

    cadd_score_runner = get_cadd_score_runner(args.db, out_dir)
    tped_runner = get_dump_tped_runner(args.db, out_dir)
    p = Pool(args.cores)
    for gene in genes:
        print "Dumping CADD scores of %s." % gene
        p.apply_async(genewise_cadd_scores, (args.db, gene, out_dir))
        print "Dumping TPED of variants in genes of %s." % gene
        p.apply_async(dump_tped, (args.db, gene, out_dir))
    p.close()
    p.join()
    print "Finished."
    # for gene in genes:
    #     genewise_cadd_scores(args.db, gene, out_dir)
    #     dump_tped(args.db, gene, out_dir)

