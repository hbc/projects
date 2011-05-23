"""Evaluate mixed populations within sequencing runs.
"""
import csv
import collections

def compare_files(call_file, expected_file, offset=0, ignore_space=False):
    """Perform comparison between a mixed call file and expected percents.
    """
    call_info = {}
    for k, v in _read_call_file(call_file, ignore_space):
        call_info[k] = v
    return compare_calls(call_info, _read_call_file(expected_file, ignore_space),
                         offset)

def compare_calls(calls, expected, offset=0):
    """Provide statistics on overlaps between calls and expected values.

    calls is a dictionary with keys of chromosomes/spaces and positions,
    and values that are dictionary of base and frequency.
    expected is a list (or generator) with chromosome/space positions and
    dictionaries of bases and frequencies.

    Returns a dictionaries with counts for:

    - single_pos: single base position with correct call
    - single_neg: single base position with wrong single call
    - single_neg_multi: single base position with multiple call
    - single_neg_multi_nomatch: single base position with multiple calls, unmatching
    - multi_pos: multiple base positions with correct call
    - multi_neg: multiple base position with wrong multi call
    - multi_neg_single: multiple base position with single call
    - multi_neg_single_nomatch: multiple base position with single call, unmatching
    """
    offset = _manage_offset_dict(offset)
    counts = collections.defaultdict(int)
    for (space, pos), ebases in expected:
        if isinstance(offset, dict):
            cur_offset = offset[pos]
            if cur_offset is None:
                continue
        else:
            cur_offset = offset
        cbases = calls.get((space, pos+cur_offset),
                           collections.defaultdict(lambda: None))
        #print space, pos, ebases, cbases
        if _is_single(ebases):
            if _is_single(cbases):
                outcome = _compare_single(ebases, cbases)
            else:
                ext = "_nomatch" if not _single_multi_match(ebases, cbases) else ""
                outcome = "single_neg_multi%s" % ext
        else:
            if _is_single(cbases):
                #print space, pos, ebases, cbases
                ext = "_nomatch" if not _single_multi_match(ebases, cbases) else ""
                outcome = "multi_neg_single%s" % ext
            else:
                outcome = _compare_multi(ebases, cbases)
        counts[outcome] += 1
    return dict(counts)

def _manage_offset_dict(offset):
    """Handle complicated offsets passed in as dictionaries.
    """
    if isinstance(offset, dict):
        point_offset = {}
        for k, val in offset.iteritems():
            start, end = k.split("-")
            for v in range(int(start), int(end)+1):
                point_offset[v] = val
        return point_offset
    else:
        return int(offset)

def _present_bases(bases):
    return [b for (b, v) in bases.iteritems() if v is not None]

def _has_base(base, bases):
    return base in _present_bases(bases)

def _compare_single(expect, call):
    if _has_base(_present_bases(expect)[0], call):
        return "single_pos"
    else:
        return "single_neg"

def _compare_multi(expect, call):
    if _present_bases(expect) == _present_bases(call):
        return "multi_pos"
    else:
        return "multi_neg"

def _is_single(bases):
    """Determine if this is a single or multiple base position.
    """
    return len(_present_bases(bases)) == 1

def _single_multi_match(one, two):
    return len(set(_present_bases(one)) & (set(_present_bases(two)))) > 0

# ## Input of read files

def _read_call_line(parts, ignore_space):
    base_order = ["A", "C", "G", "T"]
    def _parse_base(val):
        if val:
            val = float(val)
            if val <= 1:
                val *= 100.0
            if int(val) == 0:
                val = None
        else:
            val = None
        return val
    if len(parts) == 5:
        assert ignore_space
        parts = [""] + parts
    elif len(parts) == 6:
        if ignore_space:
            parts[0] = ""
    else:
        raise ValueError(parts)
    bases = dict(zip(base_order, [_parse_base(v) for v in parts[2:]]))
    return (parts[0], int(parts[1])), bases

def _read_call_file(in_file, ignore_space=False):
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        reader.next() # header
        for parts in reader:
            yield _read_call_line(parts, ignore_space)
