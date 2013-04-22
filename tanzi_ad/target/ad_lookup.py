import biomartpy as bm
from pandas import merge
import sys
from pandas.io.parsers import read_table
import sh
import os

"""
mart_name = 'ensembl'
dataset = 'hsapiens_gene_ensembl'
filters = {"entrezgene": ['54708']}

attributes = ["entrezgene", "ensembl_gene_id", "description", "hgnc_symbol"]
df = bm.make_lookup(mart_name, dataset, attributes=attributes, filters=filters)
"""

def lookup_human(filters):
    mart_name = 'ensembl'
    dataset = 'hsapiens_gene_ensembl'
    attributes = ["entrezgene", "ensembl_gene_id", "hgnc_symbol",
                  "chromosome_name", "start_position", "end_position",
                  "strand"]
    attributes = ["hgnc_symbol", "chromosome_name", "start_position",
                  "end_position", "strand"]
    df = bm.make_lookup(mart_name, dataset, attributes=attributes, filters=filters)
    return df

def write_subset(df, identifier):
    columns = ['chromosome_name', 'start_position', 'end_position']
    out_name = identifier + ".bed"
    subbed = df[df["class"] == identifier]
    subbed = subbed.ix[:, columns]
    subbed = subbed.sort(columns=columns, axis=0)
    subbed.to_csv(out_name, sep="\t", header=False, index=False)
    return out_name

def tabix(filename):
    sh.bgzip(filename)
    tab = sh.tabix.bake(p="bed")
    tab(filename + ".gz")
    return filename + ".gz"

def separate_regions(df):
    regions = df[["region" in x for x in df["symbol"]]]
    return regions

def lookup_regions(df):
    pass

def gemini_annotate(filename, db):
    base, _ = os.path.splitext(filename)
    tabixed = tabix(filename)
    annotate = sh.gemini.annotate.bake(f=tabixed, c=base, t="boolean")
    annotate(db)

if __name__ == "__main__":
    # only keep the main annotated bits not the patches
    keep_chromosome = map(str, range(1, 23) + ["MT", "X", "Y"])
    in_file = sys.argv[1]
    df = read_table(in_file, index_col=False, parse_dates=False, sep=",")
    # this file is organized with the non-early/late genes annotated at the top
    # so drop the later rows if they duplicate
    df = df.drop_duplicates(cols="symbol", take_last=False)

    # XXX IMPLEMENT THIS PART
    regions = separate_regions(df)
    regions = lookup_regions(regions)

    filters = {"hgnc_symbol": list(df["symbol"])}
    bmed = lookup_human(filters)
    bmed = bmed[[x in keep_chromosome for x in bmed["chromosome_name"]]]
    merged = merge(df, bmed, left_on="symbol", right_index=True)

    out_files = []
    out_files.append(write_subset(merged, "early_onset"))
    out_files.append(write_subset(merged, "late_onset"))
    out_files.append(write_subset(merged, "alzgene"))
    # write all
    columns = ['chromosome_name', 'start_position', 'end_position']
    master = merged.ix[:, columns]
    master = master.sort(columns=columns, axis=0)
    master.to_csv("ad_candidate.bed", sep="\t", header=False, index=False)
    out_files.append("ad_candidate.bed")
    map(tabix, out_files)
    #    map(gemini_annotate, out_files, [db] * len(out_files))

# gemini annotate -f early_onset.bed.gz -c early_onset -t boolean 1k.db
# gemini annotate -f late_onset.bed.gz -c late_onset -t boolean 1k.db
# gemini annotate -f alzgene.bed.gz -c alzgene -t boolean 1k.db
# gemini annotate -f ad_candidate.bed.gz -c ad_candidate -t boolean 1k.db
