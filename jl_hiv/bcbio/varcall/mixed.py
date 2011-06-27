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

def call_expected_iter(call_file, expected_file, offset=0, ignore_space=True):
    """Generate call and expected information for each position.
    """
    Call = collections.namedtuple('Call', ['space', 'pos', 'expected', 'called'])
    offset_pos = _pos_with_offset(offset)
    calls = {}
    for k, v in _read_call_file(call_file, ignore_space):
        calls[k] = v
    for (space, pos), ebases in _read_call_file(expected_file, ignore_space):
        opos = offset_pos(pos)
        cbases = calls.get((space, opos),
                           collections.defaultdict(lambda: None))
        yield Call(space, pos, ebases, cbases)

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
    Train = collections.namedtuple('Train', ['space', 'pos', 'base', 'freq', 'outcome'])
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    vrn_values = collections.defaultdict(list)
    train_positions = []
    offset_pos = _pos_with_offset(offset)
    for (space, pos), ebases in expected:
        opos = offset_pos(pos)
        if opos is None: continue
        cbases = calls.get((space, opos),
                           collections.defaultdict(lambda: None))
        percent_target = min(v for v in ebases.values() if v is not None)
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
        if outcome == "correct":
            call_pct, expect_pct, expect_base = _call_vrn_percent(ebases, cbases)
            train_positions.append(Train(space, opos, expect_base,
                                         expect_pct, "true_positive"))
            vrn_values[percent_target].append(call_pct)
        else:
            (_, expect_pct, _) = _call_vrn_percent(ebases, cbases)
            call_pct, call_base = sorted((v, b) for (b, v) in cbases.iteritems()
                                         if v is not None)[0]
            if call_pct < expect_pct:
                report_pct = call_pct
                wrong_type = "false_positive"
            else:
                report_pct = expect_pct
                wrong_type = "false_negative"
            train_positions.append(Train(space, opos, call_base, report_pct, wrong_type))
    return _convert_to_dict(counts), dict(vrn_values), train_positions

def _convert_to_dict(counts):
    out = {}
    for percent, vals in counts.iteritems():
        out[percent] = {}
        for k, v in vals.iteritems():
            out[percent][k] = v
    return out

def _call_vrn_percent(expect, call):
    """Retrieve the percentage variation of the call.
    """
    # first get the base we want to query -- lowest percent in the expect
    vals = [(v, b) for (b, v) in expect.iteritems() if v is not None]
    vals.sort()
    base = vals[0][1]
    return call[base], expect[base], base

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

def _pos_with_offset(offset=0):
    offset = _manage_offset_dict(offset)
    def _calc(pos):
        if isinstance(offset, dict):
            cur_offset = offset[pos]
        else:
            cur_offset = offset
        if cur_offset is not None:
            return pos+cur_offset
        else:
            return None
    return _calc

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
    def _parse_base(val, all_vals):
        total = sum(float(v) for v in all_vals if v)
        if val:
            val = float(val)
            if total < 1.1:
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
    bases = dict(zip(base_order, [_parse_base(v, parts[2:]) for v in parts[2:]]))
    return (parts[0], int(parts[1])), bases

def _read_call_file(in_file, ignore_space=False):
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        reader.next() # header
        for parts in reader:
            yield _read_call_line(parts, ignore_space)
