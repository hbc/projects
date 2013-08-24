# given a directory name of BED, count the overlaps between all pairwise
# comparisons of peaks

import pybedtools
from argparse import ArgumentParser
import glob
import os
from itertools import combinations
from pandas import DataFrame
import yaml


def count_overlaps(a, b):
    ba = pybedtools.BedTool(a)
    bb = pybedtools.BedTool(b)
    return ba.intersect(bb).count()


def extract_barcode(in_file):
    return in_file.split("_")[3]


def jaccard_index(a, b):
    ba = pybedtools.BedTool(a).sort()
    bb = pybedtools.BedTool(b).sort()
    return ba.jaccard(bb)['jaccard']


def is_same_virus(pair, config):
    virus = [config["barcodes"][x][1] for x in pair]
    return virus[0] == virus[1]


def main(args):
    directory = args.directory
    in_files = glob.glob(os.path.join(directory, "*.bed"))
    pairs = list(combinations(in_files, 2))
    overlaps = [count_overlaps(x[0], x[1]) for x in pairs]
    jaccard = [jaccard_index(x[0], x[1]) for x in pairs]
    barcode_pairs = [map(extract_barcode, x) for x in pairs]
    column_names = [x[0] + "_vs_" + x[1] for x in barcode_pairs]

    with open(args.config) as in_handle:
        config = yaml.load(in_handle)

    virus_same = [is_same_virus(x, config) for x in barcode_pairs]

    df = DataFrame({'overlap': overlaps, 'jaccard': jaccard,
                    'virus_same': virus_same},
                   index=column_names)
    grouped = df.groupby(['virus_same']).mean()
    df.to_csv("overlaps.txt", sep="\t")
    grouped.to_csv("same_virus.txt", sep="\t")

#    print "\t".join(column_names)
#    print "\t".join(map(str, overlaps))
#    print "\t".join(map(str, jaccard))


if __name__ == "__main__":
    parser = ArgumentParser("count overlaps between all BED files in a "
                            "directory")
    parser.add_argument("directory", help="directory of BED files.")
    parser.add_argument("config", help="config file to lookup barcodes")
    parser.add_argument("--barcode", help="code output table by barcode.")
    args = parser.parse_args()
    main(args)
