#!/usr/bin/env python
"""Setup file and install script for deploying the stem cell discovery engine.
"""
from setuptools import setup, find_packages

setup(name = "bcbio-scde-deploy",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "Automated Galaxy and BII installation",
      license = "MIT",
      url = "http://discovery.hsci.harvard.edu/",
      namespace_packages = ["bcbio"],
      packages = find_packages(),
      scripts = [],
      package_data = {
          'config' : ['*.yaml'],
          },
      install_requires = [
          "PyYAML >= 3.09",
          "fabric >= 1.0.1",
          "boto >= 2.0b4",
      ])
