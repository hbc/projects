#!/usr/bin/env python
"""Examine call overlaps between different barcoded samples.

Usage:
 call_overlap.py <list of files>
"""
import os
import sys
import csv

def main(*call_files):
    file_calls = {}
    for fname in call_files:
        file_calls[fname] = read_call_file(fname)
    for fname1 in call_files:
        other_files = list(call_files)
        other_files.remove(fname1)
        compare_files(fname1, other_files, file_calls)

def compare_files(fname1, other_files, file_calls):
    print "**", os.path.basename(fname1).split("-", 1)[0]
    c1 = file_calls[fname1]
    other_sets = [file_calls[f] for f in other_files]
    shared = len(c1.intersection(*other_sets))
    unique = len(c1.difference(*other_sets))
    partial_shared = len(c1) - shared - unique
    _print_w_percent("Shared all", shared, len(c1))
    _print_w_percent("Shared partial", partial_shared, len(c1))
    _print_w_percent("Unique", unique, len(c1))

def _print_w_percent(name, num, total):
    print "| % 15s | %s | %.1f%% |" % (name, num, float(num) / total * 100.0)

def read_call_file(fname):
    calls = []
    with open(fname) as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for line in reader:
            calls.append(tuple(line[:3]))
    return set(calls)

if __name__ == "__main__":
    main(*sys.argv[1:])
