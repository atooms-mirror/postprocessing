#!/usr/bin/env python

import os
import sys
import numpy
import atooms.trajectory as trj


class PrettyList(list):
    def __str__(self):
        return ' '.join([str(x) for x in self])


class TrajectoryNeighbors(trj.TrajectoryXYZ):

    """Neighbors trajectory."""

    def __init__(self, filename, mode='r', offset=1):
        super(TrajectoryNeighbors, self).__init__(filename, mode=mode)
        self._offset = offset # neighbors produced by voronoi are indexed from 1
        self.fields = ['neighbors']
        # This is necessary to format integer numpy array correctly
        self._fields_float = False


def get_neighbors(fileinp, fileout, args, fmt=None):
    if args.neigh_file is None:
        write_neighbors(fileinp, args.rcut, fileout, fmt)
        tn = trj.TrajectoryNeighbors(fileout)
    else:
        # TODO: is we ever get to make this a clean factory, test can be avoided
        if args.neigh_voronoi:
            from atooms.voronoi import TrajectoryVoronoi
            tn = TrajectoryVoronoi(args.neigh_file)
        else:
            tn = trj.TrajectoryNeighbors(args.neigh_file)

    if args.neigh_limit is not None:
        def limited(system, max_neighbors):
            for i in range(len(system.neighbors)):
                numn = len(system.neighbors[i])
                system.neighbors[i] = s.neighbors[i][0: min(max_neighbors, numn)]
            return system
        tn.register_callback(limited, args.neigh_limit)
    return tn

def compute_neighbors(system, rcut):
    import neighbors_wrap
    npart = len(system.particle)
    nmax = 300
    nn = numpy.zeros(npart, dtype=numpy.int32)
    neigh = numpy.zeros((npart, nmax), dtype=numpy.int32, order='F')
    rcut = numpy.asarray(rcut)
    pos = system.dump('pos').transpose() # copy is not needed, order not needed.
    ids = [p.id for p in system.particle]
    box = system.cell.side
    neighbors_wrap.neighbors(box, pos, ids, rcut, nn, neigh)
    for j, p in enumerate(system.particle):
        p.neighbors = neigh[j, 0:nn[j]]
    return system

def write_neighbors(fileinp, rcut, fileout='/dev/stdout', fmt=None):
    import copy
    import neighbors_wrap
    with trj.Trajectory(fileinp, fmt=fmt) as t, \
         TrajectoryNeighbors(fileout, 'w') as tout:
        tout.timestep = t.timestep
        for i, s in enumerate(t):
            s = compute_neighbors(s, rcut)
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

    
