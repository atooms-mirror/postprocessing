#!/usr/bin/env python

import os
import glob

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(name='postprocessing',
      version='0.1',
      description='Post processing tools for molecular simulations',
      author='Daniele Coslovich',
      author_email='daniele.coslovich@univ-montp2.fr',
      packages=['postprocessing'],
      scripts=glob.glob(os.path.join('bin', '*.py'))
     )
