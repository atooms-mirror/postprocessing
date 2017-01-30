#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Common neighbor analysis.

Bonds are identified by a signature (i,j,k), where 
i: particle species (1 by default)
j: number of common neighbors
k: number of of bonds between common neighbors

For instance, icosahedral bonds are 1,5,5.
"""

import os
import sys
import argparse
import numpy
from collections import defaultdict

from atooms.trajectory import Trajectory, TrajectoryNeighbors
from atooms.utils import add_first_last_skip, fractional_slice
from pyutils.histogram import Histogram
from postprocessing.neighbors import get_neighbors

parser = argparse.ArgumentParser()
parser = add_first_last_skip(parser, what=['first', 'last'])
parser.add_argument('-n',              dest='neigh_file', help='neighbors file')
parser.add_argument('-N', '--neighbor',dest='neigh', type=str, default='', help='flags for neigh.x command')
parser.add_argument('-V', '--neighbor-voronoi',dest='neigh_voronoi', action='store_true', help='neigh_file is of Voronoi type')
parser.add_argument('-M', '--neighbor-max',dest='neigh_limit', type=int, default=None, help='take up to *limit* neighbors (assuming they are ordered)')
parser.add_argument(      '--rcut', dest='rcut', help='cutoff radii as comma separated string, ex r11,r12,r22')
parser.add_argument('-o', '--output',dest='output', action='store_true', help='write to file')
parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
parser.add_argument(nargs='+', dest='files',type=str, help='input files')
args = parser.parse_args()

if len(args.tag) > 0:
    args.tag = '_' + args.tag

# Unpack cut off radii as numpy array
if args.neigh_file is None:
    if args.rcut is None:
        raise ValueError('provide cut off radii')

    rc = args.rcut.split(',')
    # This is the number of species given the number of independent cut off radii
    nsp = (-1 + int(numpy.ceil((1+8*len(rc))**0.5))) / 2
    args.rcut = numpy.ndarray((nsp, nsp))
    i = 0
    for isp in range(nsp):
        for jsp in range(isp,nsp):
            args.rcut[isp, jsp] = float(rc[i])
            i+=1

for finp in args.files:
    t = Trajectory(finp)
    tn, desc = get_neighbors(finp, args, os.path.basename(sys.argv[0]))

