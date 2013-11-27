#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name = "bcbio-hiv",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "",
      description = "",
      license = "MIT",
      url = "http://bcbio.wordpress.com",
      packages = find_packages(),
      scripts = [],
      package_data = {
          'config' : ['*.yaml']},
      install_requires = [
          "biopython >= 1.57",
          "pysam >= 0.4.1",
          "khmer>=0.7.1",
          "PyYAML >= 3.09",
          "numpy >= 1.5.1",
          "bcbio-nextgen >= 0.4a"])
