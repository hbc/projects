#!/usr/bin/env python
"""Pick targets for secondary verification from a set of top differentially expressed.

Usage:
  pick_verification_targets.py <main_target_csv> <summarized_target_csv>
"""
import sys
import os
import csv

def main(target_csv, summary_csv):
    out_csv = apply("{}-verification{}".format, os.path.splitext(target_csv))
    with open(target_csv) as in_handle:
        targets, target_name = read_targets(csv.reader(in_handle, dialect="excel-tab"))
    with open(summary_csv) as in_handle:
        targets = add_combined_info(targets, csv.reader(in_handle, dialect="excel-tab"),
                                    target_name)
    targets = filter(lambda x: x["combined_rank"] > 0, targets.values())
    targets = sorted(targets, key=lambda x: x["combined_rank"])
    with open(out_csv, "w") as out_handle:
        write_output(targets, csv.writer(out_handle, dialect="excel-tab"))

def write_output(targets, writer):
    writer.writerow(["shrnaid", "rank", "direction", "total.targets", "has.background",
                     "has.low", "biotype", "acc", "genesymbol", "description"])
    for t in targets:
        writer.writerow([t["shrnaid"],
                         t["combined_rank"],
                         "up" if t["is_up"] else "down",
                         t["total_targets"],
                         t["has_control"],
                         t["has_low"],
                         t["biotype"],
                         t["acc"],
                         t["genesymbol"],
                         t["description"]])

def add_combined_info(targets, reader, target_name):
    target_low = target_name.replace("_H", "_L")
    header = reader.next()
    for i, (acc, _, exps) in enumerate(x[:3] for x in reader):
        exps = exps.split(";")
        if target_name in exps:
            targets[acc]["combined_rank"] = i
            targets[acc]["has_low"] = target_low in exps
            targets[acc]["total_targets"] = len(exps)
    return targets

def read_targets(reader):
    header = reader.next()
    target_name = header[1]
    targets = {}
    for i, (shrna, tval, cval, _, _, acc, biotype, genesymbol, description) in \
            enumerate(x[:9] for x in reader):
        tval = float(tval)
        cval = float(cval)
        targets[acc] = {"rank": i,
                        "is_up": tval > cval,
                        "has_control": cval > 0,
                        "combined_rank": -1,
                        "shrnaid": shrna,
                        "acc": acc,
                        "biotype": biotype,
                        "genesymbol": genesymbol,
                        "description": description}
    return targets, target_name

if __name__ == "__main__":
    main(*sys.argv[1:])
