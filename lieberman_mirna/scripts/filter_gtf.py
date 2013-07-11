import os
from pandas.io.parsers import read_csv
from pandas import notnull
import pprint
import sys

DATA_DIR = "/n/home05/kirchner/hsph/projects/lieberman_mirna/results/focused_peak_calling"
IMPACT_TARGET_FILE = os.path.join(DATA_DIR, "impact_targets.txt")
GTF_FILE = "/n/home05/kirchner/hsph/biodata/genomes/Hsapiens/hg19/rnaseq/ref-transcripts.gtf"

def extract_attribute(l, attribute):
    return l.split(attribute)[1].split(';')[0].strip().replace('"', '')

def main():
    df = read_csv(IMPACT_TARGET_FILE, sep="\t")
    target_genes = list(df[notnull(df['id'])]['id'])

    with open(GTF_FILE) as in_handle:
        for line in in_handle:
            if extract_attribute(line, 'gene_id') in target_genes:
                sys.stdout.write(line)

if __name__ == "__main__":
    main()
