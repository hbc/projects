"""Multiprocessing read entry points.
"""
from bcbio import utils

from bcbio.hbc.linker import tasks
from bcbio.pipeline import lane

@utils.map_wrap
def hbc_process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def trim_with_aligner(*args):
    return tasks.trim_with_aligner(*args)
