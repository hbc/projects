from argparse import ArgumentParser

def capture_recapture(marked, captured):
    M = len(marked)
    C = len(captured)
    R = len(marked.intersection(captured))
    return float(M * C) / float(R)


if __name__ == "__main__":
    parser = ArgumentParser("Calculate population size from two samples using the "
                            "capture-recapture method")
    parser.add_argument("marked", help="IDs of initial marked population.")
    parser.add_argument("captured", help="IDs of captured population.")
    args = parser.parse_args()
    with open(args.marked) as marked_handle:
        marked = set([line.strip() for line in marked_handle])
    with open(args.captured) as captured_handle:
        captured = set([line.strip() for line in captured_handle])

    print capture_recapture(marked, captured)
