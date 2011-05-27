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

    Returns a dictionary keyed with the minimal expected percent in the expected
    calls. Each dictionary has three count values:

    - correct: Same call in the the test and expected.
    - partial: Partially correct. In the case of single (100%) regions, the
      test call had a single call plus additional calls. In the case of multiple
      expected calls, the test call had only one of the expected.
    - wrong: The test call is totally wrong compared to expected.

    So both partial and wrong are, well, wrong, but partial lets you know you are
    in the ballpark of making a correct call.
    """
    offset = _manage_offset_dict(offset)
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for (space, pos), ebases in expected:
        if isinstance(offset, dict):
            cur_offset = offset[pos]
            if cur_offset is None:
                continue
        else:
            cur_offset = offset
        cbases = calls.get((space, pos+cur_offset),
                           collections.defaultdict(lambda: None))
        percent_target = min(v for v in ebases.values() if v is not None)
        if percent_target == 99.5:
            print ebases.values()
            raise NotImplementedError
        if _is_single(ebases):
            if _is_single(cbases):
                outcome = _compare_single(ebases, cbases)
            else:
                outcome = "partial" if _single_multi_match(ebases, cbases) else "wrong"
        else:
            if _is_single(cbases):
                outcome = "partial" if _single_multi_match(ebases, cbases) else "wrong"
            else:
                outcome = _compare_multi(ebases, cbases)
        counts[percent_target][outcome] += 1
    return _convert_to_dict(counts)

def _convert_to_dict(counts):
    out = {}
    for percent, vals in counts.iteritems():
        out[percent] = {}
        for k, v in vals.iteritems():
            out[percent][k] = v
    return out

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
        return "correct"
    else:
        return "wrong"

def _compare_multi(expect, call):
    if _present_bases(expect) == _present_bases(call):
        return "correct"
    else:
        return "wrong"

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
            if val == 0.0:
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
