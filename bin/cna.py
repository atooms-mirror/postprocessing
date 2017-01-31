#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Common neighbor analysis (CNA).

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

from atooms.trajectory import Trajectory
from atooms.utils import add_first_last_skip, fractional_slice
from pyutils.histogram import Histogram
from postprocessing.neighbors import get_neighbors

def cna(particle, neighbors):
    # Add neighbor list to particles as sets
    for p, n in zip(particle, neighbors):
        p.neighbors = set(n)

    # For all i-j pairs find mutual CNA index
    data = []
    for i in range(len(particle)):
        ni = particle[i].neighbors
        for j in ni:
            if j<=i: continue
            nj = particle[j].neighbors
            # Mutual neighbors of i-j pair
            common = ni & nj
            # Count bonds between mutual neighbors
            bonds = 0
            for k in common:
                for m in common:
                    if m<=k: continue
                    if m in particle[k].neighbors:
                        bonds+=1
            # Accumulate CNA index
            #print 'cna (%s,%s): %d_%d_%d' % (i,j,1,len(common),bonds)
            data.append('%d_%d_%d' % (1,len(common),bonds))
    return data

def main(args):
    if len(args.tag) > 0:
        args.tag = '_' + args.tag

    # Unpack cut off radii as numpy array
    # TODO: refactor, this is in common with neigh.py
    if args.neigh_file is None:
        desc = 'rcut:' + args.rcut
        rc = args.rcut.split(',')
        # This is the number of species given the number of independent cut off radii
        nsp = (-1 + int(numpy.ceil((1+8*len(rc))**0.5))) / 2
        args.rcut = numpy.ndarray((nsp, nsp))
        i = 0
        for isp in range(nsp):
            for jsp in range(isp,nsp):
                args.rcut[isp, jsp] = float(rc[i])
                args.rcut[jsp, isp] = float(rc[i])
                i+=1
    else:
        desc = 'neighbors_file:' + args.neigh_file

    # Handle multiple signatures: make a list of them
    if args.signature is not None:
        args.signature = args.signature.split(',')

    for finp in args.files:
        # We need to open the trajectory here anyway, we should we do that again in get_neighbors()?
        t = Trajectory(finp, fmt=args.fmt)
        # We write neighbors to a tmp file
        import tempfile
        fout = tempfile.mkstemp()[1]
        tn = get_neighbors(finp, fout, args, fmt=args.fmt)

        # If required, we put CNA data in a separate directory
        if args.dirout is not None:
            dirout = os.path.dirname(args.dirout + '/' + finp)
            from pyutils.utils import mkdir
            mkdir(dirout)
            fbase = args.dirout + '/' + finp
        else:
            fbase = finp

        # Fraction of selected CNA bonds (signature argument)
        if args.signature is not None:
            fh = dict()
            for sign in args.signature:
                fout = fbase + '.cna%s.fraction-%s' % (args.tag, sign)
                fh[sign] = open(fout, 'w', buffering=0)
                fh[sign].write('# columns: step, fraction of CNA bond %s; %s\n' % (sign, desc))

        # Loop over samples
        hist = defaultdict(int)
        for i, s in enumerate(t):
            if t.steps[i] in tn.steps:
                ii = tn.steps.index(t.steps[i])
                data = cna(t[i].particle, tn[ii].neighbors)
                # Dump signature
                if args.signature is not None:
                    for sign in args.signature:
                        x = len([d for d in data if d == sign]) / float(len(data))
                        fh[sign].write('%d %s\n' % (t.steps[i], x))
                # Fill histogram
                for d in data:
                    hist[d]+=1

        # Write histogram
        with open(fbase + '.cna%s.hist' % args.tag, 'w') as fhhist:
            norm = sum(hist.values())
            fhhist.write('# columns: CNA bond, average; %s\n' % desc)
            for d in sorted(hist, key=hist.get, reverse=True):
                fhhist.write('%s %g\n' % (d, hist[d] / float(norm)))

        if args.signature is not None:
            for fhi in fh.values():
                fhi.close()

        if args.neigh_file is None:
            os.remove(tn.filename)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = add_first_last_skip(parser, what=['first', 'last'])
    parser.add_argument('-n',              dest='neigh_file', help='neighbors file')
    parser.add_argument('-N', '--neighbor',dest='neigh', type=str, default='', help='flags for neigh.x command')
    parser.add_argument('-V', '--neighbor-voronoi',dest='neigh_voronoi', action='store_true', help='neigh_file is of Voronoi type')
    parser.add_argument('-M', '--neighbor-max',dest='neigh_limit', type=int, default=None, help='take up to *limit* neighbors (assuming they are ordered)')
    parser.add_argument(      '--rcut', dest='rcut', help='cutoff radii as comma separated string, ex r11,r12,r22')
    parser.add_argument('-s', '--signature',dest='signature', help='signature')
    parser.add_argument('-o', '--output',dest='output', action='store_true', help='write to file')
    parser.add_argument('-d', '--dirout',dest='dirout', help='output dir')
    parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
    parser.add_argument(      '--fmt',     dest='fmt', help='input file format')
    parser.add_argument(nargs='+', dest='files',type=str, help='input files')
    args = parser.parse_args()

    main(args)
