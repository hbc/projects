#!/usr/bin/env python
"""Share BII data files with Galaxy using data libraries.

Usage:
  bii_datasets_to_galaxy.py <bii data dir>
"""
import os
import sys
import glob
import collections

from bcbio import isatab

def main(bii_dir):
    for study_dirs in bii_studies(bii_dir):
        study_datafiles(study_dirs["isatab"])

def study_datafiles(isatab_dir):
    """Retrieve data files and associated metadata for a study.
    """
    rec = isatab.parse(isatab_dir)
    assert len(rec.studies) == 1
    study = rec.studies[0]
    for assay in study.assays:
        info = _get_assay_info(assay)
        print info

def _get_assay_info(assay):
    """Retrieve all data files and samples associated with the assay.
    """
    out = collections.defaultdict(set)
    attrs = ["Raw Data File", "Derived Data File", "Sample Name"]
    for data_file, node in assay.nodes.items():
        for attr in attrs:
            for val in node.metadata.get(attr, []):
                if val:
                    out[attr].add(val)
    final = {}
    for k, v in out.iteritems():
        final[k] = list(v)
    return final

def bii_studies(bii_dir):
    """Retrieve meta data and raw data directories in BII.
    """
    bii_dirs = filter(os.path.isdir,
                      [os.path.join(bii_dir, d) for d in os.listdir(bii_dir)])
    metadata_dir = [d for d in bii_dirs if os.path.basename(d) == "meta_data"][0]
    bii_dirs.remove(metadata_dir)
    for study_dir in sorted(os.path.join(metadata_dir, d)
                            for d in os.listdir(metadata_dir)):
        study_name = os.path.basename(study_dir)
        data_dirs = filter(os.path.isdir,
                           [os.path.join(d, study_name) for d in bii_dirs])
        assert len(data_dirs) > 0, study_name
        yield dict(isatab=study_dir, data=data_dirs)

if __name__ == "__main__":
    main(*sys.argv[1:])
