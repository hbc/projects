"""Add information from a FAM or PED input file into an existing gemini database.

Whoops -- I forgot to add my PED file during load? No problem.

Usage:
   add_fam_to_gemini.py <gemini_db> <fam file>
"""

import os
import glob
import sys
import sqlite3

from gemini import ped

def main(gemini_db_file, fam_file):
    fam_info = read_fam_file(fam_file)
    for gemini_db in glob.glob(gemini_db_file):
        add_to_database(gemini_db, fam_info)

def add_to_database(db, fam_info):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "SELECT sample_id, name FROM samples"
    for row_id, name in conn.execute(sql):
        if name in fam_info:
            f = fam_info[name]
            sql = ("UPDATE samples SET family_id = ?, paternal_id = ?, maternal_id = ?, "
                   "sex = ?, phenotype = ? WHERE sample_id = ?")
            cursor.execute(sql, (f["family_id"], f["paternal_id"], f["maternal_id"], f["sex"], f["affected"], row_id))
    conn.commit()
    cursor.close()

def read_fam_file(in_file):
    out = {}
    header = ped.get_ped_fields(in_file)
    with open(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                field = line.split(None, 8)[:8]
                if len(field) > 1:
                    info = dict(zip(header, field))
                    out[info["name"]] = info
    return out


if __name__ == "__main__":
    main(*sys.argv[1:])
