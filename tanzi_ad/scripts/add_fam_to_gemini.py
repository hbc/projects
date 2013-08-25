"""Add information from a FAM or PED input file into an existing gemini database.

Whoops -- I forgot to add my PED file during load? No problem.

Usage:
   add_fam_to_gemini.py <gemini_db> <fam file>
"""

import os
import sys
import sqlite3

from gemini.ped import pedformat

def main(gemini_db, fam_file):
    add_to_database(gemini_db, read_fam_file(fam_file))

def add_to_database(db, fam_info):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    sql = "SELECT sample_id, name FROM samples"
    for row_id, name in conn.execute(sql):
        if name in fam_info:
            f = fam_info[name]
            sql = ("UPDATE samples SET family_id = ?, paternal_id = ?, maternal_id = ?, "
                   "sex = ?, phenotype = ? WHERE sample_id = ?")
            cursor.execute(sql, (f.family, f.paternal, f.maternal, f.sex, f.phenotype, row_id))
    conn.commit()
    cursor.close()

def read_fam_file(in_file):
    out = {}
    with open(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                field = line.split(None, 7)[:7]
                if len(field) > 1:
                    ped = pedformat(field)
                    out[ped.name] = ped
    return out


if __name__ == "__main__":
    main(*sys.argv[1:])
