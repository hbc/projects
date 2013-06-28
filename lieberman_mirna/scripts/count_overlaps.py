#
import sys
from collections import Counter, defaultdict

def main(in_file):
    genes = defaultdict(list)
    with open(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            read_id = line.split("=")[1].split(";")[0]
            gene_name = line.split("=")[2].split(";")[0].replace('"', "")
            genes[gene_name].append(read_id)
    for gene, reads in genes.items():
        print "\t".join([gene, str(len(set(reads)))])

if __name__ == "__main__":
    main(sys.argv[1])
