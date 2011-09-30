"""Install script for linker trimming work.
"""
from setuptools import setup, find_packages

setup(name = "bcbio-hbc-linker",
      version = "0.1",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "Dual side linker trimming",
      license = "MIT",
      url = "https://github.com/chapmanb/bcbb/tree/master/nextgen",
      namespace_packages = ["bcbio", "bcbio.hbc"],
      packages = find_packages(),
      scripts = ['scripts/hsph_trim_by_align.py'],
      install_requires = [
          "bcbio-nextgen >= 0.3a",
      ])
