import os
import subprocess
from argparse import ArgumentParser

def gemini_region(chrom, start, end, out_file):
    if file_exists(out_file):
        return out_file
    region = "{chrom}:[start}-{end}".format(**locals())
    with file_transaction(out_file) as tx_out_file:
        query_cmd = ('gemini query --format tped --region {region} '
                    '-q "select * from variants" > {tx_out_file}')
        subprocess.check_call(query_cmd.format(**locals()), shell=True)
    return out_file

if __name__ == "__main__":
    description = ("Script to dump regions from AD_Imputation_Burden_Gene_"
                   "co-ordinates.txt to TPED format with regions named by "
                   "gene.")
    parser = ArgumentParser(description=description)
    parser.add_argument("--out-dir", default=".", help="Optional output directory.")
    parser.add_argument("burden_coordinates",
                        help=("Burden coordinates file, tab delimited with the format "
                              "name chr start end status"))
    args = parser.parse_args()

    with open(args.burden_coordinates) as in_handle:
        # skip header
        in_handle.next()
        for line in in_handle:
            name, chrom, start, end, _ = line.split("\t")
            out_file = os.path.join(args.out_dir, name) + ".tped"
            print "Processing {chrom}, {start}, {end}.".format(**locals())
            #gemini_region(chrom, start, end, out_file)
