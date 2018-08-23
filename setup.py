#!/usr/bin/env python

import os
import glob

# We use numpy distutils to compile and wrap f90 code via f2py
from numpy.distutils.core import setup, Extension

with open('atooms/postprocessing/_version.py') as f:
    exec(f.read())

args = dict(name='atooms-pp',
            version=__version__,
            description='Post-processing tools for molecular simulations',
            author='Daniele Coslovich',
            author_email='daniele.coslovich@umontpellier.fr',
            url='http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/',
            packages=['atooms/postprocessing'],
            scripts=glob.glob(os.path.join('bin', '*.py')),
            install_requires=['atooms>=1', 'numpy', 'argh'],
            ext_modules=[Extension('postprocessing.realspace_wrap',
                                   sources=['atooms/postprocessing/realspace.f90']),
                         Extension('postprocessing.fourierspace_wrap', 
                                   sources=['atooms/postprocessing/fourierspace.f90'])],
            license='GPLv3',
            classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Development Status :: 5 - Production/Stable',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 2',
                'Topic :: Scientific/Engineering :: Physics',
            ]
)

setup(**args)
