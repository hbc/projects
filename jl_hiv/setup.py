#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name = "bcbio-hiv",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "",
      description = "",
      license = "MIT",
      url = "http://bcbio.wordpress.com",
      namespace_packages = ["bcbio"],
      packages = find_packages(),
      scripts = [],
      package_data = {
          'config' : ['*.yaml']},
      install_requires = [
          "biopython >= 1.57",
          "PyYAML >= 3.09"])
