#!/usr/bin/env python
import argparse
from intermine.webservice import Service

def query(ids):
    service = Service("http://targetmine.nibio.go.jp/targetmine")
    query = service.new_query("Protein")
    query.add_view(
        "primaryIdentifier", "primaryAccession", "name", "length",
        "compounds.compound.casRegistryNumber", "compounds.compound.name",
        "compounds.compound.compoundGroup.name"
    )
    test_id = ids[0]
    query.add_constraint("Protein", "IN", ",".join(ids))
    return query.rows()

def _get_ids(input):
    ids = input.next()
    ids = ids.strip().split(" ")
    return ids

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type = argparse.FileType('r'), default = '-')
    args = parser.parse_args()
    ids = _get_ids(args.input)

    query_rows = query(ids)
    for row in query_rows:
        print row["primaryIdentifier"], row["primaryAccession"], row["name"], row["length"],\
             row["compounds.compound.casRegistryNumber"], row["compounds.compound.name"],\
            row["compounds.compound.compoundGroup.name"]
