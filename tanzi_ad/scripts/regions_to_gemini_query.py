"""
Convert a csv of chrom, start, ref, alt to a GEMINI query string
"""
from argparse import ArgumentParser
from gemini import GeminiQuery


def parse_csv_line(line):
    chrom, pos, ref, alt = line.split(",")
    return chrom, pos, ref, alt

def csv_to_query_string(csv_handle, db):
    query_list = []
    count = 0
    for line in csv_handle:
        chrom, pos, ref, alt = line.strip().split(",")
        query_list.append("(chrom='{chrom}' and start={pos})".format(**locals()))
        count += 1
        if count % 100 == 0:
            gq = GeminiQuery(db)
            regions = " or ".join(query_list)
            query = "select * from variants where %s" % regions
            gq.run(query)
            if count == 100:
                print gq.header
            for row in gq:
                print row
            query_list = []

if __name__ == "__main__":
    description = ("Convert csv of chrom, start, ref, alt to a GEMINI "
                   "query string.")
    parser = ArgumentParser(description=description)
    parser.add_argument("CSV", help="CSV file of variants")
    parser.add_argument("db", help="GEMINI database.")
    args = parser.parse_args()

    with open(args.CSV) as csv_handle:
        csv_to_query_string(csv_handle, args.db)
