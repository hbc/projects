"""
takes a fasta file and outputs the appropriate adaptor sequence
arguments for cutadapt
"""
import fileinput
import sys

def main():
    for line in fileinput.input():
        # ignore fasta headers
        if line[0] == ">":
            continue
        if line.strip():
            print "-b " + "\"" + line.strip() + "\"",

if __name__ == "__main__":
    main()
