#!/usr/bin/env python

import os
import glob

# We use numpy distutils to compile and wrap f90 code via f2py
from numpy.distutils.core import setup, Extension

args = dict(name='postprocessing',
            version='0.1',
            description='Post processing tools for molecular simulations',
            author='Daniele Coslovich',
            author_email='daniele.coslovich@umontpellier.fr',
            url='http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/',
            packages=['postprocessing'],
            scripts=glob.glob(os.path.join('bin', '*.py')),
            install_requires=['atooms>=1', 'numpy', 'scipy', 'argh'],
            ext_modules=[Extension('postprocessing.neighbors_wrap', 
                                   sources=['postprocessing/neighbors.f90'])]
)

setup(**args)
