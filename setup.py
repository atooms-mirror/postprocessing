#!/usr/bin/env python

import os
import glob
from setuptools import setup, find_packages

setup(name='postprocessing',
      version=0.1,
      description='Post processing tools for molecular simulations',
      author='Daniele Coslovich',
      author_email='daniele.coslovich@univ-montp2.fr',
      packages=find_packages(exclude=('tests', 'docs')),
      install_requires=['argh'],
      scripts=glob.glob(os.path.join('bin', '*.py'))
     )
