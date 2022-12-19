# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Linked cells to compute neighbors efficiently."""

from collections import defaultdict
import numpy

# TODO: define iterator over cells

class LinkedCells(object):

    def __init__(self, rcut):
        self.rcut = rcut
        self.neighbors = []
        self._is_adjusted = False

    def adjust(self, box, newton, periodic=None):
        self.box = box
        self.hbox = box / 2
        self.n_cell = (box / self.rcut).astype(int)
        self.box_cell = self.box / self.n_cell
        if periodic is None:
            periodic = numpy.array([True] * len(box))
        self._map(newton, periodic)

    def _map(self, newton, periodic):
        def _pbc(t, N):
            """Apply PBCs to cell only along periodic directions"""
            for i in range(len(t)):
                if not periodic[i]:
                    continue
                if t[i] >= N[i]:
                    t[i] -= N[i]
                elif t[i] < 0:
                    t[i] += N[i]
            return t

        def _outside(t, N):
            """Check whether subcell is outside cell"""
            for i in range(len(t)):
                if t[i] >= N[i] or t[i] < 0:
                    return True
            return False

        def _map_3d_newton(n_cell):
            neigh_cell = {}
            for ix in range(n_cell[0]):
                for iy in range(n_cell[1]):
                    for iz in range(n_cell[2]):
                        # This is a map of neighbouring cells obeying III law
                        neigh_cell[(ix, iy, iz)] = \
                            [(ix+1, iy, iz), (ix+1, iy+1, iz), (ix, iy+1, iz),
                             (ix-1, iy+1, iz), (ix+1, iy, iz-1), (ix+1, iy+1, iz-1),
                             (ix, iy+1, iz-1), (ix-1, iy+1, iz-1), (ix+1, iy, iz+1),
                             (ix+1, iy+1, iz+1), (ix, iy+1, iz+1), (ix-1, iy+1, iz+1),
                             (ix, iy, iz+1)]
            return neigh_cell

        def _map_2d_newton(n_cell):
            neigh_cell = {}
            for ix in range(n_cell[0]):
                for iy in range(n_cell[1]):
                    # This is a map of neighbouring cells obeying III law
                    neigh_cell[(ix, iy)] = \
                        [(ix+1, iy), (ix+1, iy+1), (ix, iy+1), (ix-1, iy+1)]
            return neigh_cell

        def _map_3d_nonewton(n_cell):
            neigh_cell = {}
            for ix in range(n_cell[0]):
                for iy in range(n_cell[1]):
                    for iz in range(n_cell[2]):
                        neigh_cell[(ix, iy, iz)] = []
                        for deltax in [-1, 0, 1]:
                            for deltay in [-1, 0, 1]:
                                for deltaz in [-1, 0, 1]:
                                    if deltax == deltay == deltaz == 0:
                                        continue
                                    neigh_cell[(ix, iy, iz)].append((ix+deltax, iy+deltay, iz+deltaz))
            return neigh_cell

        def _map_2d_nonewton(n_cell):
            neigh_cell = {}
            for ix in range(n_cell[0]):
                for iy in range(n_cell[1]):
                    neigh_cell[(ix, iy)] = []
                    for deltax in [-1, 0, 1]:
                        for deltay in [-1, 0, 1]:
                            if deltax == deltay == 0:
                                continue
                            neigh_cell[(ix, iy)].append((ix+deltax, iy+deltay))
            return neigh_cell

        if len(self.n_cell) == 3:
            if newton:
                self._neigh_cell = _map_3d_newton(self.n_cell)
            else:
                self._neigh_cell = _map_3d_nonewton(self.n_cell)
        elif len(self.n_cell) == 2:
            if newton:
                self._neigh_cell = _map_2d_newton(self.n_cell)
            else:
                self._neigh_cell = _map_2d_nonewton(self.n_cell)
        else:
            raise ValueError('linked cells not supported for dimensions not in {2,3}')

        # Apply PBC
        for key in self._neigh_cell:
            for idx in range(len(self._neigh_cell[key])):
                folded = _pbc(list(self._neigh_cell[key][idx]), self.n_cell)
                self._neigh_cell[key][idx] = tuple(folded)

        # Remove subcells that are out of bounds
        # (this is only applied to non periodic directions)
        new_neigh_cell = {}
        for key in self._neigh_cell:
            new_neigh_cell[key] = []
            for subcell in self._neigh_cell[key]:
                if not _outside(subcell, self.n_cell):
                    new_neigh_cell[key].append(subcell)
        self._neigh_cell = new_neigh_cell

    def _index(self, pos):
        x = ((pos + self.hbox) / self.box_cell)
        return x.astype(numpy.int)

    def on_border(self, pos):
        index = list(self._index(pos))
        found = False
        for i in range(len(index)):
            if index[i] == self.n_cell[i] - 1 or index[i] == 0:
                found = True
                break
        return found

    def compute(self, box, pos, other=None, as_array=False, newton=True, periodic=None):
        if not self._is_adjusted:
            self.adjust(box, newton=other is None and newton, periodic=periodic)
            self._is_adjusted = True
        # We only need positions here but how can we be sure that
        # this is the same set of particles we use when retrieving
        # the neighbours? We should keep a reference.
        # from atooms.postprocessing.realspace_wrap import compute
        # index = numpy.ndarray(pos.shape, dtype=numpy.int)
        # compute.bin_in_cell(pos, self.hbox, self.box_cell, index)

        self.neighbors = []
        self.number_of_neighbors = []
        index = self._index(pos)
        if other is None:
            index_other = index
        else:
            index_other = self._index(other)

        particle_in_cell = defaultdict(list)
        for ipart, icell in enumerate(index_other):
            particle_in_cell[tuple(icell)].append(ipart)

        for ipart in range(pos.shape[0]):
            icell = tuple(index[ipart])
            # Initialize an empty list
            neighbors = []
            # Start with particles in the cell of particle ipart
            if other is None:
                neighbors += [_ for _ in particle_in_cell[icell] if _ > ipart]
            else:
                neighbors += particle_in_cell[icell]
                # try:
                #     neighbors.remove(ipart)
                # except:
                #     pass
            # Loop over neighbors cells and add neighbors
            for jcell in self._neigh_cell[icell]:
                neighbors += particle_in_cell[jcell]
            self.neighbors.append(neighbors)
            self.number_of_neighbors.append(len(neighbors))

        if as_array:
            npart = len(self.neighbors)
            number_of_neighbors = numpy.array(self.number_of_neighbors)
            neighbors_array = numpy.ndarray((npart, max(number_of_neighbors)), dtype=numpy.int)
            for ipart in range(len(self.neighbors)):
                neighbors_array[ipart, 0:len(self.neighbors[ipart])] = self.neighbors[ipart]
            self.neighbors = neighbors_array
            self.number_of_neighbors = number_of_neighbors
            return self.neighbors, number_of_neighbors
        else:
            return self.neighbors

# import sys
# from atooms.system import Particle
# import atooms.trajectory as trj
# t = trj.Trajectory(sys.argv[1])
# system = t[0]
# # system.particle = []
# # dr = system.cell.side[0] / 20
# # for ix in range(0,20):
# #     for iy in range(0,20):
# #         system.particle.append(Particle(position=[ix*dr+dr/2 - system.cell.side[0]/2, iy*dr+dr/2 - system.cell.side[0]/2, 0]))

# trj.decorators.change_species(system, 'F')
# print('----')
# pos = system.dump('pos', order='C')
# ids = system.dump('spe')
# print('----')
# lc = _LinkedCells(rcut=2.0)
# nn, num = lc.compute(system.cell.side, pos, as_array=True)
# print(nn[0][0:num[0]])
# print('#', lc.box_cell / (system.cell.side/2))
# print('----')

# for j in nn[0][0:num[0]]:
#     dr = system.particle[0].distance(system.particle[j], system.cell)
#     print(numpy.sum(dr**2)**0.5, lc.box_cell[0] * 2)

# # from atooms.postprocessing.realspace_wrap import compute
# # rcut = numpy.array([[1.5, 1.5], [1.5, 1.5]])
# # rcut[:] = lc.box_cell[0] * 2
# # max_neighbors = 300
# # number_of_neighbors = numpy.ndarray(pos.shape[0], dtype=numpy.int32)
# # neighbors = numpy.ndarray((pos.shape[0], max_neighbors), dtype=numpy.int32, order='F')
# # #print(pos.shape)
# # pos = numpy.array(pos, order='F')
# # compute.neighbors_list('C',system.cell.side,pos.transpose(),ids,rcut,number_of_neighbors,neighbors)
# # #neighbors_0 = numpy.ndarray(max_neighbors, dtype=numpy.int32, order='F')
# # #nn = numpy.array(0, dtype=numpy.int32)
# # # compute.neighbors('C',system.cell.side,pos,pos[0],ids,rcut,nn,neighbors_0)
# # # print(neighbors[0][0:number_of_neighbors[0]])


# # print(str(system.particle[0].position / (system.cell.side/2))[1:-1])
# # print()
# # print()
# # for nnn in [sorted(neighbors[0]), sorted(lc.neighbors[0])]:
# #     for j in nnn:
# #         print(str(system.particle[j].position / (system.cell.side/2))[1:-1])
# #     print()
# #     print()


# # js = []
# # for j in range(pos.shape[0]):
# #     dr = system.particle[0].distance(system.particle[j], system.cell)
# #     if j> 0 and numpy.sum(dr**2)**0.5 < lc.box_cell[0] * 2:
# #         js.append(j)

# # # print(len(sorted(lc.neighbors[0])))
# # # print(len(sorted(neighbors[0][0:number_of_neighbors[0]])))
# # # print(len(sorted(js)))
# # # print('----')

# # # for nnn in [sorted(lc.neighbors[0]),
# # #             sorted(neighbors[0][0:number_of_neighbors[0]])]:
# # #     for j in nnn:
# # #         print(str(system.particle[j].position)[1:-1])
# # #     print()
# # #     print()

# # # for nnn in [sorted(js),
# # #             sorted(neighbors[0][0:number_of_neighbors[0]])]:
# # #     for j in nnn:
# # #         dr = system.particle[0].distance(system.particle[j], system.cell)
# # #         print(j, numpy.sum(dr**2)**0.5)
# # #     print


# # # # for x in neighbors:
# # # #     print(x)
# # # for i in range(pos.shape[0]):
# # #     x = list(sorted(lc.neighbors[i]))
# # #     y = neighbors[i][0:number_of_neighbors[i]]
# # #     print(x)
# # #     print(y)
# # #     assert numpy.all(x == y)
# # #     print()
