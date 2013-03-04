import MySQLdb as mdb
import sys
import logging as logger
import unittest
import os
import argparse


def _load_file(in_file, header=False):
    rsids = []
    locs = []
    with open(in_file) as in_handle:
        # skip the header
        if header:
            in_handle.readline()
        for line in in_handle:
            if "rs" in line:
                rsids.append(line)
            else:
                locs.append(line)
    return rsids, locs


def _convert_locations(locs):
    converted = []
    for loc in locs:
        spl = [x.strip() for x in loc.split("-")]
        converted.append(str(spl[0]) + ":" + str(spl[1]) + "-" + str(spl[1]))

    return converted


def _clean_rsids(data):
    def _format_rsid(x):
        return '"' + x.strip() + '"'
    stripped = [_format_rsid(x) for x in data]
    return ", ".join(stripped)


def _lookup_rsids(con, table, rsids):
    cur = con.cursor()
    cur.execute('SELECT * from %s WHERE name in (%s)' % (table, rsids))
    return cur


def _lookup_chunk(con, table, rsids, chunk_size=100):
    left = rsids[chunk_size:]
    cleaned = _clean_rsids(rsids[0:chunk_size])
    cur = _lookup_rsids(con, table, cleaned)
    for snp in cur:
        chrom = snp[1].replace("chr", "")
        print chrom + ":" + str(snp[2]) + "-" + str(snp[2]) + "\t" +  snp[4]
    return left


def _connect_db(user, password, database, url='localhost'):
    try:
        con = mdb.connect(url, user, password, database)
    except mdb.Error, e:
        logger.error("Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit(1)
    return con


def main(in_file, user, password, table, db):
    rsids, locs = _load_file(in_file, header=False)
    locs = _convert_locations(locs)
    #rsids = _clean_rsids(rsids)
    con = _connect_db(user, password, db)
    while rsids:
        rsids = _lookup_chunk(con, table, rsids)
    for loc in locs:
        print loc.strip()


class TestConverter(unittest.TestCase):

    TEST_FILE = os.path.join(os.path.dirname(__file__), "test_snps.txt")

    def test_load_file(self):
        rsids, locs = _load_file(self.TEST_FILE, header=True)
        self.assertTrue(all([len(rsids) == 5, len(locs) == 4]))

    def test_convert_locations(self):
        locs = ["1-731181\n", "1-751456\n"]
        correct = 'chr1\t731181\n'
        converted = _convert_locations(locs)
        self.assertEquals(correct, converted[0])

    def test_clean_rsids(self):
        rsids = ["rs12565286\n", "rs11804171\n"]
        cleaned = _clean_rsids(rsids)
        correct = '"rs12565286", "rs11804171"'
        self.assertEquals(correct, cleaned)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Return locations of rsids "
                                     "by looking them up in a mysql db.")
    parser.add_argument('--user')
    parser.add_argument('--password')
    parser.add_argument('--db')
    parser.add_argument('--table')
    parser.add_argument('in_file')

    args = parser.parse_args()
    main(args.in_file, args.user, args.password, args.db, args.table)
