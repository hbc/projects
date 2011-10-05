"""Multiprocessing read entry points.
"""
from bcbio import utils

from bcbio.hbc.linker import tasks

@utils.map_wrap
def trim_with_aligner(*args):
    return tasks.trim_with_aligner(*args)
