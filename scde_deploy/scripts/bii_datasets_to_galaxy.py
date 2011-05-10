#!/usr/bin/env python
"""Share BII data files with Galaxy using data libraries.

Usage:
  bii_datasets_to_galaxy.py <YAML config file> <Study config file>
"""
import os
import sys
import glob
import operator
import collections
import subprocess

import yaml

from bcbio import isatab
from bcbio.galaxy.api import GalaxyApiAccess

def main(config_file, study_config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with open(study_config_file) as in_handle:
        study_config = yaml.load(in_handle)
    bii_dir = os.path.join(config["base_install"], config["bii_dirname"],
                           config["bii_data_dirname"])
    groups = [("Study Assay Technology Type", "Study Assay Measurement Type"),
              ("organism",),
              ("organism part",),
              ("Study Title",)]
    prepped_files = organize_datafiles(bii_dir, groups, study_config)
    add_to_galaxy_datalibs(prepped_files, config)

# ## Synchronize BII files with Galaxy data libraries

def add_to_galaxy_datalibs(prepped_files, config):
    """Add the organized files to synchronized Galaxy data libraries.

    3 actions needed:
      - create data library for each top level item
      - create folders and sub-folders to subsequent levels
      - add links to data in final folders
    """
    galaxy_api = GalaxyApiAccess(config["galaxy_url"], config["galaxy_apikey"])
    for key, vals in prepped_files.iteritems():
        dl_name = "SCDE: %s" % (", ".join([k for k in key if k]))
        _add_data_library(galaxy_api, dl_name, vals)

def _add_data_library(galaxy_api, name, info):
    print name
    data_library = galaxy_api.get_datalibrary_id(name)
    cur_items = galaxy_api.library_contents(data_library)
    root = [f for f in cur_items
            if f["type"] == "folder" and f["name"] == "/"][0]
    for folder, info in _add_folders(info, root, data_library,
                                     cur_items, galaxy_api):
        _add_data_links(data_library, folder, info, cur_items, galaxy_api)

def _add_folders(info, parent_folder, data_library, cur_items,
                 galaxy_api):
    """Recursively add folders to get at study items.
    """
    out = []
    for fname, items in info.iteritems():
        folder = _add_or_get_folder(data_library, parent_folder, fname[0],
                                    cur_items, galaxy_api)
        if isinstance(items, dict):
            out.extend(_add_folders(items, folder, data_library, cur_items, galaxy_api))
        else:
            out.append((folder, items))
    return out

def _add_or_get_folder(data_library, base_folder, new_name, cur_items, galaxy_api):
    """Retrieve or add the new folder
    """
    folder_path = os.path.join(base_folder["name"], new_name)
    existing = [f for f in cur_items
                if f["type"] == "folder" and f["name"] == folder_path]
    if len(existing) == 0:
        folder = galaxy_api.create_folder(data_library, base_folder["id"], new_name)
        return folder[0]
    else:
        assert len(existing) == 1
        return existing[0]

def _add_data_links(data_library, folder, items, cur_items, galaxy_api):
    """Provide links to the raw and processed BII data from within Galaxy.
    """
    def _get_ext(data_file):
        return os.path.splitext(data_file["name"])[-1].lower()
    def _gunzip(fname):
        out_fname = os.path.splitext(fname)[0]
        if not os.path.exists(out_fname):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
        return out_fname
    act_exts = {".gz": _gunzip}
    want_exts = [".tiff", ".cel", ".txt", "", ".bed"]
    org_builds = {"Homo sapiens (Human)" : "hg19",
                  "Homo sapiens": "hg19",
                  "Mus musculus (Mouse)": "mm9",
                  "Mus musculus musculus": "mm9",
                  "Rattus norvegicus (Rat)": "rn4"}
    data_types = collections.defaultdict(list)
    for x in items:
        data_types[x["type"].split()[0]].append(x)
    for data_type, data_items in data_types.iteritems():
        data_folder = _add_or_get_folder(data_library, folder, data_type,
                                         cur_items, galaxy_api)
        for data_file in data_items:
            ext = _get_ext(data_file)
            if act_exts.has_key(ext):
                data_file["name"] = act_exts[ext](data_file["name"])
                ext = _get_ext(data_file)
            if ext in want_exts or len(ext) > 10:
                full_path = os.path.join(data_folder["name"],
                                         os.path.basename(data_file["name"]))
                existing = [f for f in cur_items
                            if f["name"] == full_path and f["type"] == "file"]
                if len(existing) == 0:
                    galaxy_api.upload_from_filesystem(data_library, data_folder["id"],
                                                      data_file["name"],
                                                      org_builds[data_file["organism"][0]])
            else:
                print data_file["name"]

# ## Organize data files into high level structure for Galaxy deployment

def organize_datafiles(bii_dir, groups, study_config):
    all_files = []
    for study_dirs in bii_studies(bii_dir):
        for data_file in study_datafiles(study_dirs["isatab"],
                                         study_dirs["data"], study_config):
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
        name = tuple([item.get(g, [""])[0] for g in groups])
        out[name].append(item)
    return dict(out)

# ## Parse ISA-Tab metadata, retrieving files of interest

def study_datafiles(isatab_dir, data_dirs, study_config):
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
    attrs = [("organism", "taxID", "Organism"),
             ("strain", "developmental stage", "health status", "phenotype"),
             ("organism part", "Organism Part", "cell type", "primary tumor site"),
             ]
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
            keys = node.metadata.keys()
            vals = [node.metadata[k][0] for k in keys]
            raise ValueError("Missing %s %s %s" % (attr_group, keys, vals))
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
