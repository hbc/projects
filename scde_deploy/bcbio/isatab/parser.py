"""Parse ISA-Tab structured metadata describing experimental data.

Works with ISA-Tab (http://isatab.sourceforge.net), which provides a structured
format for describing experimental metdata.
"""
from __future__ import with_statement

import os
import csv
import glob
import collections

def parse(isatab_ref):
    """Entry point to parse an ISA-Tab directory.

    isatab_ref can point to a directory of ISA-Tab data, in which case we
    search for the investigator file, or be a reference to the high level
    investigation file.
    """
    if os.path.isdir(isatab_ref):
        fnames = glob.glob(os.path.join(isatab_ref, "i_*.txt"))
        assert len(fnames) == 1
        isatab_ref = fnames[0]
    assert os.path.exists(isatab_ref), "Did not find investigation file: %s" % isatab_ref
    i_parser = InvestigationParser()
    with open(isatab_ref) as in_handle:
        rec = i_parser.parse(in_handle)
    s_parser = StudyAssayParser(isatab_ref)
    rec = s_parser.parse(rec)
    return rec

class InvestigationParser:
    """Parse top level investigation files into ISATabRecord objects.
    """
    def __init__(self):
        self._sections = {
            "ONTOLOGY SOURCE REFERENCE": "ontology_refs",
            "INVESTIGATION": "metadata",
            "INVESTIGATION PUBLICATIONS": "publications",
            "INVESTIGATION CONTACTS": "contacts",
            "STUDY DESIGN DESCRIPTORS": "design_descriptors",
            "STUDY PUBLICATIONS": "publications",
            "STUDY FACTORS": "factors",
            "STUDY ASSAYS" : "assays",
            "STUDY PROTOCOLS" : "protocols",
            "STUDY CONTACTS": "contacts"}
        self._nolist = ["metadata"]

    def parse(self, in_handle):
        line_iter = self._line_iter(in_handle)
        # parse top level investigation details
        rec = ISATabRecord()
        rec, _ = self._parse_region(rec, line_iter)
        # parse study information
        while 1:
            study = ISATabStudyRecord()
            study, had_info = self._parse_region(study, line_iter)
            if had_info:
                rec.studies.append(study)
            else:
                break
        return rec

    def _parse_region(self, rec, line_iter):
        """Parse a section of an ISA-Tab, assigning information to a supplied record.
        """
        had_info = False
        keyvals, section = self._parse_keyvals(line_iter)
        if keyvals:
            rec.metadata = keyvals[0]
        while section and section[0] != "STUDY":
            had_info = True
            keyvals, next_section = self._parse_keyvals(line_iter)
            attr_name = self._sections[section[0]]
            if attr_name in self._nolist:
                try:
                    keyvals = keyvals[0]
                except IndexError:
                    keyvals = {}
            setattr(rec, attr_name, keyvals)
            section = next_section
        return rec, had_info

    def _line_iter(self, in_handle):
        """Read tab delimited file, handling ISA-Tab special case headers.
        """
        reader = csv.reader(in_handle, dialect="excel-tab")
        for line in reader:
            if len(line) > 0 and line[0]:
                # check for section headers; all uppercase and a single value
                if line[0].upper() == line[0] and "".join(line[1:]) == "":
                    line = [line[0]]
                yield line

    def _parse_keyvals(self, line_iter):
        """Generate dictionary from key/value pairs.
        """
        out = None
        line = None
        for line in line_iter:
            if len(line) == 1:
                break
            else:
                # setup output dictionaries, trimming off blank columns
                if out is None:
                    while not line[-1]:
                        line = line[:-1]
                    out = [{} for _ in line[1:]]
                # add blank values if the line is stripped
                while len(line) < len(out) + 1:
                    line.append("")
                for i in range(len(out)):
                    out[i][line[0]] = line[i+1]
                line = None
        return out, line

class StudyAssayParser:
    """Parse row oriented metadata associated with study and assay samples.

    This currently does not attempt to be complete, but rather to extract the
    most useful information (in my biased opinion) and represent it simply
    in the record objects.

    This is coded generally, so can be expanded to more cases. It is biased
    towards microarray and next-gen sequencing data.
    """
    def __init__(self, base_file):
        self._dir = os.path.dirname(base_file)
        self._col_quals = ("Parameter Value", "Performer", "Date", "Unit",
                           "Term Accession Number", "Term Source REF")
        self._col_types = {"attribute": ("Characteristics", "Factor Type",
                                         "Comment", "Label", "Material Type"),
                           "node" : ("Sample Name", "Source Name", "Image File",
                                     "Raw Data File", "Derived Data File"),
                           "node_assay" : ("Extract Name", "Labeled Extract Name",
                                           "Assay Name", "Data Transformat Name",
                                           "Normalization Name"),
                           "processing": ("Protocol REF")}
        self._synonyms = {"Array Data File" : "Raw Data File",
                          "Derived Array Data File" : "Derived Data File",
                          "Hybridization Assay Name": "Assay Name",
                          "Derived Array Data Matrix File": "Raw Data File",
                          "Raw Spectral Data File": "Raw Data File",
                          "Derived Spectral Data File": "Derived Data File"}

    def parse(self, rec):
        """Retrieve row data from files associated with the ISATabRecord.
        """
        for study in rec.studies:
            source_data = self._parse_study(study.metadata["Study File Name"],
                                            "Sample Name")
            for assay in study.assays:
                assay_data = self._parse_study(assay["Study Assay File Name"],
                                               "Raw Data File")
        return rec

    def _parse_study(self, fname, node_type):
        """Parse study or assay row oriented file around the supplied base node.
        """
        nodes = {}
        with open(os.path.join(self._dir, fname)) as in_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            header = self._swap_synonyms(reader.next())
            hgroups = self._collapse_header(header)
            htypes = self._characterize_header(header, hgroups)
            name_index = header.index(node_type)
            for line in reader:
                name = line[name_index]
                try:
                    node = nodes[name]
                except KeyError:
                    node = NodeRecord(name, node_type)
                attrs = self._line_keyvals(line, header, hgroups, htypes)
                print attrs
                for i, oi in enumerate(g[0] for g in hgroups):
                    print htypes[i], header[oi], line[oi]
                print
                break
        return nodes.values()

    def _line_keyvals(self, line, header, hgroups, htypes):
        out = collections.defaultdict(list)
        out = self._line_by_type(line, header, hgroups, htypes, out, "attribute")
        out = self._line_by_type(line, header, hgroups, htypes, out, "processing",
                                 self._extract_protocol_info)
        out = self._line_by_type(line, header, hgroups, htypes, out, "node")
        return dict(out)

    def _line_by_type(self, line, header, hgroups, htypes, out, want_type,
                      collapse_quals_fn = None):
        """Parse out key value pairs for line information based on a group of values.
        """
        for index, htype in ((i, t) for i, t in enumerate(htypes) if t == want_type):
            col = hgroups[index][0]
            key = header[col]
            if key.find("[") >= 0:
                key = key.replace("]", "").split("[")[-1]
            if collapse_quals_fn:
                val = collapse_quals_fn(line, header, hgroups[index])
            else:
                val = line[col]
            out[key].append(val)
        return out

    def _extract_protocol_info(self, line, header, indexes):
        print indexes
        for i in indexes:
            print header[i], line[i]
        print

    def _characterize_header(self, header, hgroups):
        """Characterize header groups into different data types.
        """
        out = []
        for h in [header[g[0]] for g in hgroups]:
            this_ctype = None
            for ctype, names in self._col_types.iteritems():
                if h.startswith(names):
                    this_ctype = ctype
                    break
            out.append(this_ctype)
        return out

    def _collapse_header(self, header):
        """Combine header columns into related groups.
        """
        out = []
        for i, h in enumerate(header):
            if h.startswith(self._col_quals):
                out[-1].append(i)
            else:
                out.append([i])
        return out

    def _swap_synonyms(self, header):
        return [self._synonyms.get(h, h) for h in header]

class ISATabRecord:
    """Represent ISA-Tab metadata in structured format.

    High level key/value data.
      - metadata -- dictionary
      - ontology_refs -- list of dictionaries
      - contacts -- list of dictionaries
      - publications -- list of dictionaries

    Sub-elements:
      - studies: List of ISATabStudyRecord objects.
    """
    def __init__(self):
        self.metadata = {}
        self.ontology_refs = []
        self.publications = []
        self.contacts = []
        self.studies = []

class ISATabStudyRecord:
    """Represent a study within an ISA-Tab record.
    """
    def __init__(self):
        self.metadata = {}
        self.design_descriptors = []
        self.publications = []
        self.factors = []
        self.assays = []
        self.protocols = []
        self.contacts = []

class NodeRecord:
    """Represent a data node within an ISA-Tab Study/Assay file.
    """
    def __init__(self, name="", ntype=""):
        self.ntype = ntype
        self.name = name
        self.metadata = {}
