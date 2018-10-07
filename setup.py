#!/usr/bin/env python

import os
import glob

# We use numpy distutils to compile and wrap f90 code via f2py
import setuptools
from numpy.distutils.core import setup, Extension


# Get the long description from README.md and try to convert it to
# reST. Adapted from https://bons.ai/blog/markdown-for-pypi
try:
    from pypandoc import convert
    readme = convert('README.md', 'rst')
except (ImportError, OSError):
    try:
        readme = open('README.md', 'r').read()
    except:
        readme = ''

with open('atooms/postprocessing/_version.py') as f:
    exec(f.read())

args = dict(name='atooms-pp',
            version=__version__,
            description='Post-processing tools for particle simulations',
            long_description=readme,
            author='Daniele Coslovich',
            author_email='daniele.coslovich@umontpellier.fr',
            url='http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/',
            packages=['atooms', 'atooms/postprocessing'],
            scripts=glob.glob(os.path.join('bin', '*.py')),
            install_requires=['atooms>=1', 'numpy', 'argh', 'tqdm'],
            ext_modules=[Extension('atooms.postprocessing.realspace_wrap',
                                   sources=['atooms/postprocessing/realspace.f90']),
                         Extension('atooms.postprocessing.fourierspace_wrap', 
                                   sources=['atooms/postprocessing/fourierspace.f90'])],
            license='GPLv3',
            setup_requires = ['numpy'],
            classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Development Status :: 5 - Production/Stable',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 2',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.4',
                'Topic :: Scientific/Engineering :: Physics',
            ]
)

setup(**args)
