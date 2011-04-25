"""Tests for parsing and extracting information from ISA-Tab formatted metadata.
"""
import os
import unittest

from bcbio import isatab

class IsatabTest(unittest.TestCase):
    def setUp(self):
        self._dir = os.path.join(os.path.dirname(__file__), "isatab")

    def test_basic_parsing(self):
        """Test general parsing of an example ISA directory.
        """
        work_dir = os.path.join(self._dir, "BII-I-1")
        rec = isatab.parse(work_dir)
        assert rec.metadata["Investigation Identifier"] == "BII-I-1"
        assert len(rec.ontology_refs) == 6
        assert rec.ontology_refs[2]["Term Source Name"] == "UO"
        assert len(rec.publications) == 1
        assert rec.publications[0]["Investigation Publication DOI"] == "doi:10.1186/jbiol54"

        assert len(rec.studies) == 2
        study = rec.studies[0]
        assert study.metadata["Study File Name"] == "s_BII-S-1.txt"
        assert len(study.assays) == 3
        assert study.assays[0]["Study Assay File Name"] == "a_metabolome.txt"

    def test_minimal_parsing(self):
        """Parse a minimal ISA-Tab file without some field values filled in.
        """
        work_dir = os.path.join(self._dir, "minimal")
        rec = isatab.parse(work_dir)
        assert len(rec.publications) == 0
        assert len(rec.metadata) == 0

        assert len(rec.studies) == 1
        assert len(rec.studies[0].design_descriptors) == 0
