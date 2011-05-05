#!/usr/bin/env python
"""Share BII data files with Galaxy using data libraries.

Usage:
  bii_datasets_to_galaxy.py <YAML config file>
"""
import os
import sys
import glob
import operator
import collections

import yaml

from bcbio import isatab

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    bii_dir = os.path.join(config["base_install"], config["bii_dirname"],
                           config["bii_data_dirname"])
    groups = [("Study Assay Technology Type", "Study Assay Measurement Type"),
              ("organism",),
              ("organism part",),
              ("Study Title",)]
    prepped_files = organize_datafiles(bii_dir, groups)
    add_to_galaxy_datalibs(prepped_files, len(groups))

# ## Synchronize BII files with Galaxy data libraries

def add_to_galaxy_datalibs(prepped_files, levels):
    """Add the organized files to synchronized Galaxy data libraries.
    """
    pass

# ## Organize data files into high level structure for Galaxy deployment

def organize_datafiles(bii_dir, groups):
    all_files = []
    for study_dirs in bii_studies(bii_dir):
        for data_file in study_datafiles(study_dirs["isatab"],
                                         study_dirs["data"]):
            all_files.append(data_file)
    organized_files = _organize_by_first_group(all_files, groups)
    cur = [organized_files]
    for i in range(len(groups)):
        for c in cur:
            print c.keys(), len(c)
        cur = reduce(operator.add, [c.values() for c in cur])
    for c in cur:
        print c[0]
    return organized_files

def _organize_by_first_group(items, groups):
    """Recursively organize items based on subsetting into groups.
    """
    cur_group = _group_by(items, groups[0])
    if len(groups) == 1:
        return cur_group
    else:
        out = {}
        for key, vals in cur_group.iteritems():
            int_group = _organize_by_first_group(vals, groups[1:])
            out[key] = int_group
        return out

def _group_by(items, groups):
    """Group a set of items based on the provided groups.
    """
    out = collections.defaultdict(list)
    for item in items:
        name = tuple([item[g][0] for g in groups])
        out[name].append(item)
    return dict(out)

# ## Parse ISA-Tab metadata, retrieving files of interest

def study_datafiles(isatab_dir, data_dirs):
    """Retrieve data files and associated metadata for a study.
    """
    ftypes = ["Derived Data File", "Raw Data File"]
    rec = isatab.parse(isatab_dir)
    assert len(rec.studies) == 1
    study = rec.studies[0]
    for assay in study.assays:
        study_info = _get_study_matadata(study, assay)
        info = _get_assay_info(assay)
        verbose = study.metadata["Study Identifier"] in ["SB-S-29"]
        sample_info = _get_sample_metadata(study, info["Sample Name"],
                                           verbose)
        study_info.update(sample_info)
        for ftype in ftypes:
            for fname in info.get(ftype, []):
                try:
                    out = {"name": _get_full_path(fname, ftype, data_dirs),
                           "type": ftype}
                    out.update(study_info)
                    yield out
                except ValueError:
                    print "Missing file", fname, data_dirs

def _get_full_path(name, ftype, data_dirs):
    """Retrieve the full path in BII to the data file.
    """
    ftype_dirs = {"Derived Data File": "processed_data",
                  "Raw Data File": "raw_data"}
    for data_dir in data_dirs:
        test_fname = os.path.join(data_dir, ftype_dirs[ftype], name)
        if os.path.exists(test_fname):
            return test_fname
    raise ValueError("Did not find %s in %s" % (name, data_dirs))

def _get_study_matadata(study, assay):
    attrs = ["Study Assay Technology Type", "Study Assay Technology Platform",
             "Study Assay Measurement Type", "Study Title", "Study PubMed ID",
             "Study Person Last Name", "Study Identifier"]
    out = collections.defaultdict(set)
    for item in [study.metadata, assay.metadata] + study.publications + study.contacts:
        for attr in attrs:
            val = item.get(attr, "")
            if val:
                out[attr].add(val)
    return _collectionset_to_dict(out)

def _get_sample_metadata(study, samples, verbose=False):
    attrs = [("organism part", "Organism Part", "cell type"),
             ("organism", "taxID", "Organism"),
             ("strain", "developmental stage", "health status", "phenotype")]
    out = collections.defaultdict(set)
    for s in samples:
        node = study.nodes.get(s, None)
        if node:
            out = _fill_attrs_from_node(node, attrs, out)
        # we have a mismatch in the annotated sample and the actual names
        # pick a sample from the study to use instead
        else:
            out = _fill_attrs_from_node(study.nodes.values()[0], attrs, out)
    return _collectionset_to_dict(out)

def _get_assay_info(assay):
    """Retrieve all data files and samples associated with the assay.
    """
    attrs = [("Raw Data File",), ("Derived Data File",), ("Sample Name",)]
    out = collections.defaultdict(set)
    for node in assay.nodes.values():
        out = _fill_attrs_from_node(node, attrs, out)
    return _collectionset_to_dict(out)

def _fill_attrs_from_node(node, attrs, out):
    allowed_miss = ["strain", "Derived Data File"]
    for attr_group in attrs:
        found = False
        added = False
        for attr in attr_group:
            for val in node.metadata.get(attr, []):
                found = True
                if val:
                    out[attr_group[0]].add(val)
                    added = True
            if added:
                break
        if not found and attr_group[0] not in allowed_miss:
            raise ValueError("Missing %s %s" % (attr_group, node.metadata.keys()))
    return out

def _collectionset_to_dict(out):
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
