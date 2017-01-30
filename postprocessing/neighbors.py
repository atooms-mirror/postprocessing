#!/usr/bin/env python

import os
import sys
import numpy
import atooms.trajectory as trj

class Limited(object):

    # TODO: refactor it as get_neighbors

    """Limit the number of neighbors"""

    def __new__(cls, component, max_neighbors):
        cls = type('Limited', (Limited, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, max_neighbors=None):
        self.max_neighbors = max_neighbors

    def read_sample(self, sample):
        s = super(Limited, self).read_sample(sample)
        if self.max_neighbors is not None:
            for i in range(len(s.neighbors)):
                numn = len(s.neighbors[i])
                s.neighbors[i] = s.neighbors[i][0:min(self.max_neighbors, numn)]
        return s

class PrettyList(list):
    def __str__(self):
        return ' '.join([str(x) for x in self])

class TrajectoryNeighbors(trj.TrajectoryXYZ):

    """Neighbors trajectory."""

    def __init__(self, filename, mode='r', offset=1):
        super(TrajectoryNeighbors, self).__init__(filename, mode=mode)
        self._offset = offset # neighbors produced by voronoi are indexed from 1
        self.fmt = ['neighbors']
        # This is necessary to format integer numpy array correctly
        self._fmt_float = False

def get_neighbors(f, args, tag):
    if args.neigh_file is None:
        fn = f + '.%s.neigh' % tag
        compute_neighbors(f, args.rcut, fn)
        tn = trj.TrajectoryNeighbors(fn)
    else:
        from atooms.plugins.voronoi import TrajectoryVoronoi
        # TODO: is we ever get to make this a clean factory, test can be avoided
        if args.neigh_voronoi:
            from atooms.plugins.voronoi import TrajectoryVoronoi
            tn = TrajectoryVoronoi(args.neigh_file)
        else:
            tn = trj.TrajectoryNeighbors(args.neigh_file)
    if args.neigh_limit is not None:
        tnl = Limited(tn, args.neigh_limit)
    else:
        tnl = tn
    return tnl

def compute_neighbors(fileinp, rcut, fileout='/dev/stdout'):
    import copy
    import neighbors_module
    with trj.Trajectory(fileinp) as t, \
         TrajectoryNeighbors(fileout, 'w') as tout:
        npart = len(t[0].particle)
        nmax = 300
        nn = numpy.zeros(npart, dtype=numpy.int32)
        neigh = numpy.zeros((npart, nmax), dtype=numpy.int32, order='F')
        rcut = numpy.asarray(rcut)
        for i, s in enumerate(t):
            pos = s.dump('pos').transpose() # copy is not needed, order not needed.
            ids = [p.id for p in s.particle]
            box = s.cell.side
            neighbors_module.neighbors_module.neighbors(box, pos, ids, rcut, nn, neigh)
            for j, p in enumerate(s.particle):
                p.neighbors = neigh[j, 0:nn[j]]
            tout.write_sample(s, t.steps[i])

def all_neighbors(s):
    neigh = []
    for i, pi in enumerate(s.particle):
        nn = numpy.array([j for j in range(len(s.particle)) if i != j])
        neigh.append(nn)
    return neigh

def main(t, kcut=1.0):
    for s in t:
        npart = len(s.particle)
        for i in range(npart):
            for j in range(npart):
                if i <= j: continue
                dr = s.particle[i].distance(s.particle[j], s.cell)
                if numpy.sqrt(numpy.dot(dr,dr)) < kcut:
                    print i, j

    
