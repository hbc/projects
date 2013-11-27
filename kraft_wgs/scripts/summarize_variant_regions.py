"""Summarize interesting genomic differences in gemini database of population variant calls.

Retrieve gene at position, then identify nearby variations.

For usage information:
  summarize_variant_regions.py -h
"""
import argparse
import os
import sys

from gemini import GeminiQuery

def main(args):
    gene = get_impact_gene(args)
    summarize_gene_region(args, gene)

def summarize_gene_region(args, gene):
    gq = GeminiQuery(args.geminidb)
    gq.run("SELECT chrom, start, end, ref, alt, type, "
           "num_hom_ref, num_het, num_hom_alt, "
           "aaf_1kg_all, "
           "gene, impact, impact_severity, aa_change, clinvar_sig, "
           "grc, gms_illumina, in_cse, rmsk "
           "FROM variants WHERE "
           "chrom == '{chrom}' AND gene == '{gene}' AND filter is NULL "
           "ORDER BY start"
           .format(chrom=args.chrom, gene=gene))
    for row in gq:
        if row["impact_severity"] not in ["LOW"]:
            print row
            var_depths = []
            novar_depths = []
            for i, (gt, gt_type, gt_depth) in enumerate(zip(row.gts, row.gt_types, row.gt_depths)):
                if gt_type > 0:
                    print "  ", gq.index2sample[i], gt, gt_depth
        elif row["type"] == "indel":
            print row["chrom"], row["start"], row["ref"], row["alt"]

def get_impact_gene(args):
    gq = GeminiQuery(args.geminidb)
    gq.run("SELECT gene FROM variants WHERE chrom == '{chrom}' AND start == {pos}"
           .format(chrom=args.chrom, pos=args.pos - 1))
    return gq.next()["gene"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Exploratory script for summarizing population regions of interest")
    parser.add_argument("geminidb", help="Gemini database of variations",
                        type=os.path.abspath)
    parser.add_argument("chrom")
    parser.add_argument("pos", type=int)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())







