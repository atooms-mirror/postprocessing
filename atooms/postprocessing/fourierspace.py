# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Fourier-space post processing code."""

import sys
import numpy
import math
import random
import warnings
from collections import defaultdict

from atooms.trajectory.utils import check_block_size, is_cell_variable
from .helpers import linear_grid, logx_grid, adjust_skip, setup_t_grid
from .correlation import Correlation
from .qt import self_overlap


def expo_sphere(k0, kmax, pos):

    """Returns the exponentials of the input positions for each k."""

    # Technical note: we use ellipsis, so that we can pass either a
    # single sample or multiple samples without having to add a
    # trivial extra dimension to input array
    im = numpy.complex(0.0, 1.0)
    # The integer grid must be the same as the one set in kgrid,
    # otherwise there is an offset the problem is that integer
    # negative indexing is impossible in python and rounding or
    # truncating kmax can slightly offset the grid

    # We pick up the smallest k0 to compute the integer grid
    # This leaves many unused vectors in the other directions, which
    # could be dropped using different nkmax for x, y, z
    nk_max = 1 + int(kmax / min(k0))
    expo = numpy.ndarray((len(pos), ) + pos[0].shape + (2*nk_max+1, ), numpy.complex)
    expo[..., nk_max] = numpy.complex(1.0, 0.0)
    # First fill positive k
    for j in range(pos[0].shape[-1]):
        expo[..., j, nk_max+1] = numpy.exp(im * k0[j] * pos[..., j])
        expo[..., j, nk_max-1] = expo[..., j, nk_max+1].conjugate()
        for i in range(2, nk_max):
            expo[..., j, nk_max+i] = expo[..., j, nk_max+i-1] * expo[..., j, nk_max+1]
    # Then take complex conj for negative ones
    for i in range(2, nk_max+1):
        expo[..., nk_max+i] = expo[..., nk_max+i-1] * expo[..., nk_max+1]
        expo[..., nk_max-i] = expo[..., nk_max+i].conjugate()

    return expo

def expo_sphere_safe(k0, kmax, pos):
    """
    Returns the exponentials of the input positions for each k.
    It does not use ellipsis.
    """
    im = numpy.complex(0.0, 1.0)
    ndims = pos.shape[-1]
    nk_max = 1 + int(kmax / min(k0))
    expo = numpy.ndarray(pos.shape + (2*nk_max+1, ), numpy.complex)
    expo[:, :, :, nk_max] = numpy.complex(1.0, 0.0)

    for j in range(ndims):
        expo[:, :, j, nk_max+1] = numpy.exp(im*k0[j]*pos[:, :, j])
        expo[:, :, j, nk_max-1] = expo[:, :, j, nk_max+1].conjugate()
        for i in range(2, nk_max):
            expo[:, :, j, nk_max+i] = expo[:, :, j, nk_max+i-1] * expo[:, :, j, nk_max+1]

    for i in range(2, nk_max+1):
        expo[:, :, :, nk_max+i] = expo[:, :, :, nk_max+i-1] * expo[:, :, :, nk_max+1]
        expo[:, :, :, nk_max-i] = expo[:, :, :, nk_max+i].conjugate()

    return expo

def k_norm(ik, k0):
    if isinstance(k0, list) or isinstance(k0, numpy.ndarray):
        return math.sqrt((k0[0]*ik[0])**2 + (k0[1]*ik[1])**2 + (k0[2]*ik[2])**2)
    else:
        return math.sqrt(float(ik[0]**2 + ik[1]**2 + ik[2]**2)) * k0


class FourierSpaceCorrelation(Correlation):

    """
    Base class for Fourier space correlation functions.
    
    The correlation function is computed for each of the scalar values
    k_i of the provided `kgrid`. If the latter is `None`, the grid is
    built using `ksamples` entries linearly spaced between `kmin` and
    `kmax`.

    For each sample k_i in `kgrid`, the correlation function is
    computed over at most `nk` wave-vectors (k_x, k_y, k_z) such that
    their norm (k_x^2+k_y^2+k_z^2)^{1/2} lies within `dk` of the
    prescribed value k_i.
    """
    
    def __init__(self, trajectory, grid, symbol, short_name,
                 description, phasespace, nk=8, dk=0.1, kmin=-1, kmax=10,
                 ksamples=20):
        super(FourierSpaceCorrelation, self).__init__(trajectory,
                                                      grid, symbol, short_name,
                                                      description, phasespace)
        if not self._need_update:
            return

        # Some additional variables. k0 = smallest wave vectors
        # compatible with the boundary conditions
        self.nk = nk
        self.dk = dk
        self.kmin = kmin
        self.kmax = kmax
        self.ksamples = ksamples

        # Find k grid. It will be copied over to self.grid at end
        variables = self.symbol.split('(')[1][:-1]
        variables = variables.split(',')
        if len(variables) > 1:
            self.kgrid = grid[variables.index('k')]
        else:
            self.kgrid = grid

        # Setup grid once. If cell changes we'll call it again
        self._setup()

    def _setup(self, sample=0):
        self.k0 = 2*math.pi/self.trajectory[sample].cell.side
        # If grid is not provided, setup a linear grid from kmin,kmax,ksamples data
        # TODO: This shouldnt be allowed with fluctuating cells
        # Or we should fix the smallest k to some average of smallest k per sample
        if self.kgrid is None:
            if self.kmin > 0:
                self.kgrid = linear_grid(self.kmin, self.kmax, self.ksamples)
            else:
                self.kgrid = linear_grid(min(self.k0), self.kmax, self.ksamples)
        else:
            # If the first wave-vector is negative we replace it by k0
            if self.kgrid[0] < 0.0:
                self.kgrid[0] = min(self.k0)

        # Setup the grid of wave-vectors
        self.kvec, self.kvec_centered = self._setup_grid_sphere(len(self.kgrid)*[self.dk],
                                                                self.kgrid, self.k0)

    def _setup_grid_sphere(self, dk, kgrid, k0):
        """
        Setup wave vector grid with spherical average (no symmetry),
        picking up vectors that fit into shells of width dk centered around
        the values specified in the input list kgrid.
        Returns a dictonary of lists of wavevectors, one entry for each element in the grid.
        """
        kvec = defaultdict(list)
        kvec_centered = defaultdict(list)
        # With elongated box, we choose the smallest k0 component to setup the integer grid
        # This must be consistent with expo_grid() otherwise it wont find the vectors
        kmax = kgrid[-1] + dk[-1]
        kbin_max = 1 + int(kmax / min(k0))
        # TODO: it would be more elegant to define an iterator over ix, iy, iz for sphere, hemisphere, ... unless kmax is very high it might be more efficient to operate on a 3d grid to construct the vectors
        kmax_sq = kmax**2
        for ix in range(-kbin_max, kbin_max+1):
            for iy in range(-kbin_max, kbin_max+1):
                for iz in range(-kbin_max, kbin_max+1):
                    # Slightly faster and more explicit than
                    #   ksq = sum([(x*y)**2 for x, y in zip(k0, [ix, iy, iz])])
                    ksq = ((k0[0]*ix)**2 + (k0[1]*iy)**2 + (k0[2]*iz)**2)
                    if ksq > kmax_sq:
                        continue
                    # beware: numpy.sqrt is x5 slower than math one!
                    knorm = math.sqrt(ksq)
                    # Look for a shell of vectors in which the vector could fit.
                    # This expression is general and allows arbitrary k grids
                    # However, searching for the shell like this is not fast
                    # (it costs about as much as the above)
                    for ki, dki in zip(kgrid, dk):
                        if abs(knorm - ki) < dki:
                            kvec[ki].append((ix+kbin_max, iy+kbin_max, iz+kbin_max))
                            kvec_centered[ki].append((ix, iy, iz))
                            break

        # if len(kvec.keys()) != len(kgrid):
        #     _log.info('some k points could not be found')

        return kvec, kvec_centered

    def _decimate_k(self):
        """
        Pick up a random, unique set of nk vectors out ot the avilable
        ones without exceeding maximum number of vectors in shell
        nkmax.
        """
        # Setting the seed here once so as to get the same set
        # independent of filters.
        random.seed(1)
        k_sorted = sorted(self.kvec.keys())
        k_selected = []
        for knorm in k_sorted:
            nkmax = len(self.kvec[knorm])
            k_selected.append(random.sample(list(range(nkmax)), min(self.nk, nkmax)))
        return k_sorted, k_selected

    def report(self, k_sorted, k_selected):
        s = []
        for kk, knorm in enumerate(k_sorted):
            av = 0.0
            for i in k_selected[kk]:
                av += k_norm(self.kvec_centered[knorm][i], self.k0)
            s.append("# k %g : k_av=%g (nk=%d)" % (knorm, av / len(k_selected[kk]),
                                                           len(k_selected[kk])))
            # for i in k_selected[kk]:
            #     s.append('%s' % (self.kvec_centered[knorm][i] * self.k0))
        return '\n'.join(s)

    def _actual_k_grid(self, k_sorted, k_selected):
        k_grid = []
        for kk, knorm in enumerate(k_sorted):
            av = 0.0
            for i in k_selected[kk]:
                av += k_norm(self.kvec_centered[knorm][i], self.k0)
            k_grid.append(av / len(k_selected[kk]))
        return k_grid
