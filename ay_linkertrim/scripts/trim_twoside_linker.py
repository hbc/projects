#!/usr/bin/env python
"""Trim linker sequences present on both sides of insert sequence.

Usage:
  trim_twoside_linker.py <YAML config>
"""
import os, sys, unittest

import yaml

from Bio import Seq, pairwise2

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

def linker_pair_pos(seq, linker1, linker2, n_mismatch):
    """Find start and end of internal sequence flanked by a linker pair.
    """
    # sort linkers to process the longest first
    linker1, linker2 = sorted([linker1, linker2], key=len, reverse=True)
    s1, e1 = linker_pos(seq, linker1, n_mismatch)
    if s1 is not None:
        # get the largest region on either side of the linker
        seq_s, seq_e = sorted([(0, s1), (e1, len(seq))],
                              key=lambda xs: xs[1] - xs[0], reverse=True)[0]
        remain_seq = seq[seq_s:seq_e]
        s2, e2 = linker_pos(remain_seq, linker2, n_mismatch)
        if e2 is not None:
            s2 += seq_s
            e2 += seq_s
            if s2 > e1:
                return e1, s2
            else:
                return e2, s1
    # did not find pair
    return None, None

def linker_pos(seq, linker, n_mismatch):
    """Find the position of the linker with the specified number of mismatches.

    Returns None is the linker is not found in the sequence.
    """
    gap_char = "-"
    pos = str(seq).find(linker)
    if pos > 0:
        seq_region = seq[pos:pos+len(linker)]
        linker_region = linker
    else:
        aligns = pairwise2.align.localms(str(seq), str(linker),
                5.0, -4.0, -9.0, -0.5, one_alignment_only=True,
                gap_char=gap_char)
        if len(aligns) > 0:
            seq_a, linker_a, _, start, end = aligns[0]
            linker_region = linker_a[start:end]
            seq_region = seq_a[start:end]
            pos = seq.find(seq_region.replace(gap_char, ""))
        else:
            seq_region, linker_region = ("", "")
    matches = sum((1 if s == linker_region[i] else 0) for i, s in
            enumerate(seq_region))
    gaps = seq_region.count("-")
    cur_mismatch = len(linker) - matches + gaps
    if cur_mismatch <= n_mismatch:
        assert pos >= 0
        return pos, pos + len(seq_region.replace(gap_char, ""))
    else:
        return (None, None)

# ## Testing code

class TwoSideLinkerTest(unittest.TestCase):
    """Test removal of two sided linkers from sequence data.
    """
    def test_1_linker_pos(self):
        """Find a linker, with mismatches, in a larger sequence.
        """
        assert linker_pos("TTGATCTT", "GATC", 0) == (2, 6)
        assert linker_pos("TTGATCTT", "GGTC", 0) == (None, None)
        assert linker_pos("TTGATCTT", "GGTC", 1) == (2, 6)
        assert linker_pos("NATCTT", "GATC", 1) == (1, 4)
        assert linker_pos("NATCTT", "GATC", 0) == (None, None)

    def test_2_linker_pair_pos(self):
        """Find pair of linkers in a larger sequence
        """
        self.assertEqual(linker_pair_pos("TGATCTTTTCGATCGATTTT", "GATC", "CGATCGA", 0),
                         (5, 9))
        self.assertEqual(linker_pair_pos("TNATCTTTTCGATCGATTTT", "GATC", "CGATCGA", 0),
                         (None, None))
        self.assertEqual(linker_pair_pos("TNATCTTTTCGATCGATTTT", "GATC", "CGATCGA", 1),
                         (5, 9))
        self.assertEqual(linker_pair_pos("TNATCTTTTCGANCGATTTT", "GATC", "CGATCGA", 1),
                         (5, 9))
        self.assertEqual(linker_pair_pos("TCGATCGATTTTGATCTTTT", "GATC", "CGATCGA", 0),
                         (8, 12))

if __name__ == "__main__":
    main(*sys.argv[1:])
