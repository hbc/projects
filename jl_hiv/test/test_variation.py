"""Test code for calling and evaluating variations.
"""
import os
import unittest

from bcbio.variation.mixed import compare_calls

class MixedVariationEvaluation(unittest.TestCase):
    """Evaluation of variation in mixed populations.
    """
    def _percents(self, vals):
        bases = ["G", "A", "T", "C"]
        return dict(zip(bases, vals))
    def _expected(self, pos, vals):
        return [[("", pos), self._percents(vals)]]
    def _calls(self, pos, vals):
        return {("", pos): self._percents(vals)}
    def test_call_simple(self):
        """Basic mixed variation sample calling.
        """
        vals = [100.0, None, None, None]
        vals2 = [None, 100.0, None, None]
        vals3 = [None, None, 50.0, 50.0]
        vals4 = [50.0, 50.0, None, None]
        pairs = [(vals, vals), (vals2, vals), (vals4, vals), (vals3, vals),
                 (vals3, vals3), (vals4, vals3), (vals, vals4), (vals, vals3)]
        expect = ["single_pos", "single_neg", "single_neg_multi",
                  "single_neg_multi_nomatch",
                  "multi_pos", "multi_neg", "multi_neg_single",
                  "multi_neg_single_nomatch"]
        for (one, two), e in zip(pairs, expect):
            c = compare_calls(self._calls(1, one), self._expected(1, two))
            assert c.get(e, 0) == 1, (one, two, c, e)

    def test_call_offset(self):
        """Call with an offset in the expected bases.
        """
        vals = [100.0, None, None, None]
        c = compare_calls(self._calls(1, vals), self._expected(2, vals), offset=-1)
        assert c.get("single_pos", 0) == 1
        c = compare_calls(self._calls(4, vals), self._expected(2, vals), offset=2)
        assert c.get("single_pos", 0) == 1
