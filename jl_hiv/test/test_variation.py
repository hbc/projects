"""Test code for calling and evaluating variations.
"""
import os
import unittest

from bcbio.varcall.mixed import compare_calls

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
        pairs = [(vals, vals), (vals2, vals), (vals4, vals),
                 (vals3, vals), (vals3, vals3), (vals4, vals3),
                 (vals, vals4), (vals, vals3)]
        expect = [(100.0, "correct"), (100.0, "wrong"), (100.0, "partial"),
                  (100.0, "wrong"), (50.0, "correct"), (50.0, "wrong"),
                  (50.0, "partial"), (50.0, "wrong")]
        for (one, two), (e1, e2) in zip(pairs, expect):
            c = compare_calls(self._calls(1, one), self._expected(1, two))
            assert c.get(e1, {}).get(e2, 0) == 1, (one, two, c, e1, e2)

    def test_call_offset(self):
        """Call with an offset in the expected bases.
        """
        vals = [100.0, None, None, None]
        c = compare_calls(self._calls(1, vals), self._expected(2, vals), offset=-1)
        assert c.get(100.0, {}).get("correct", 0) == 1
        c = compare_calls(self._calls(4, vals), self._expected(2, vals), offset=2)
        assert c.get(100.0, {}).get("correct", 0) == 1

    def test_call_multiple(self):
        """Test for incorrect calls with multiple samples.
        """
        vals = [95.0, 5.0, None, None]
        vals2 = [95.0, 4.0, 1.0, None]

        pairs = [(vals2, vals), (vals, vals)]
        expect = [(5.0, "wrong"), (5.0, "correct")]
        for (one, two), (e1, e2) in zip(pairs, expect):
            c = compare_calls(self._calls(1, one), self._expected(1, two))
            assert c.get(e1, {}).get(e2, 0) == 1, (one, two, c, e1, e2)
