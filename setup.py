#!/usr/bin/env python

import os
import glob

args = dict(name='postprocessing',
            version='0.1',
            description='Post processing tools for molecular simulations',
            author='Daniele Coslovich',
            author_email='daniele.coslovich@umontpellier.fr',
            url='http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/',
            packages=['postprocessing'],
            scripts=glob.glob(os.path.join('bin', '*.py')))

try:
    from setuptools import setup, find_packages
    args['package_data'] = {
        'postprocessing': ['*.so', '*.f90'],
    }
except ImportError:
    from distutils.core import setup, Extension
    args['ext_modules'] = [Extension('postprocessing', ['f90wrap_neighbors.f90', '_neighbors_module.so'])],

setup(**args)
