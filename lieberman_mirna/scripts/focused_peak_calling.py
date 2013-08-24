
# extract list of genes from impact targets file
import os
from pandas.io.parsers import read_csv

DATA_DIR = "/n/home05/kirchner/hsph/projects/lieberman_mirna/results/focused_peak_calling"
IMPACT_TARGETS = os.path.join(DATA_DIR, "impact_targets.txt")

df = read_csv(IMPACT_TARGETS, sep="\t")
target_genes = list(df['id'])
print target_genes
