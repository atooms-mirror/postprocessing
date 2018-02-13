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
from .realspace import self_overlap

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
    for j in xrange(pos[0].shape[-1]):
        expo[..., j, nk_max+1] = numpy.exp(im * k0[j] * pos[..., j])
        expo[..., j, nk_max-1] = expo[..., j, nk_max+1].conjugate()
        for i in xrange(2, nk_max):
            expo[..., j, nk_max+i] = expo[..., j, nk_max+i-1] * expo[..., j, nk_max+1]
    # Then take complex conj for negative ones
    for i in xrange(2, nk_max+1):
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

    for j in xrange(ndims):
        expo[:, :, j, nk_max+1] = numpy.exp(im*k0[j]*pos[:, :, j])
        expo[:, :, j, nk_max-1] = expo[:, :, j, nk_max+1].conjugate()
        for i in xrange(2, nk_max):
            expo[:, :, j, nk_max+i] = expo[:, :, j, nk_max+i-1] * expo[:, :, j, nk_max+1]

    for i in xrange(2, nk_max+1):
        expo[:, :, :, nk_max+i] = expo[:, :, :, nk_max+i-1] * expo[:, :, :, nk_max+1]
        expo[:, :, :, nk_max-i] = expo[:, :, :, nk_max+i].conjugate()

    return expo

def k_norm(ik, k0):
    if isinstance(k0, list) or isinstance(k0, numpy.ndarray):
        return math.sqrt((k0[0]*ik[0])**2 + (k0[1]*ik[1])**2 + (k0[2]*ik[2])**2)
    else:
        return math.sqrt(float(ik[0]**2 + ik[1]**2 + ik[2]**2)) * k0


class FourierSpaceCorrelation(Correlation):

    def __init__(self, trajectory, grid, name, short_name,
                 description, phasespace, nk=8, dk=0.1, kmin=-1, kmax=10,
                 ksamples=20):
        # grid and name variables can be lists or tuples, ex. ['k', 't'] or ['k', 'w']
        # TODO: the time grid is not used here
        super(FourierSpaceCorrelation, self).__init__(trajectory,
                                                      grid, name, short_name,
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
        if type(name) is list or type(name) is tuple:
            self.kgrid = grid[name.index('k')]
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
        for ix in xrange(-kbin_max, kbin_max+1):
            for iy in xrange(-kbin_max, kbin_max+1):
                for iz in xrange(-kbin_max, kbin_max+1):
                    ksq = sum([(x*y)**2 for x, y in zip(k0, [ix, iy, iz])])
                    if ksq > kmax**2:
                        continue
                    # beware: numpy.sqrt is x5 slower than math one!
                    knorm = math.sqrt(ksq)
                    # look for a shell in which the vector could fit
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
        ones # without exceeding maximum number of vectors in shell
        nkmax.
        """
        # Setting the seed here once so as to get the same set
        # independent of filters.
        random.seed(1)
        k_sorted = sorted(self.kvec.keys())
        k_selected = []
        for knorm in k_sorted:
            nkmax = len(self.kvec[knorm])
            k_selected.append(random.sample(range(nkmax), min(self.nk, nkmax)))
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


class SelfIntermediateScattering(FourierSpaceCorrelation):

    #TODO: xyz files are 2 slower than hdf5 where
    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=8, tsamples=60,
                 dk=0.1, kmin=1.0, kmax=10.0, ksamples=10, norigins=-1):
        FourierSpaceCorrelation.__init__(self, trajectory, [kgrid, tgrid], ('k', 't'), \
                                         'fkt.self', 'Self intermediate scattering function',
                                         'pos-unf', nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        # Before setting up the time grid, we need to check periodicity over blocks
        check_block_size(self.trajectory.steps, self.trajectory.block_size)
        if tgrid is None:
            self.grid[1] = [0.0] + logx_grid(trajectory.timestep,
                                             trajectory.total_time * 0.75, tsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid[1])
        self.skip = adjust_skip(trajectory, norigins)

        # TODO: Can this be moved up?
        self.k_sorted, self.k_selected = self._decimate_k()
        #self.log.info(self.report(self.k_sorted, self.k_selected))
        #self.log.debug(self.k_selected)

    def _compute(self):
        # Throw everything into a big numpy array
        # TODO: remove this redundancy
        self._pos = self._pos_unf
        pos = numpy.ndarray((len(self._pos), ) + self._pos[0].shape)
        for i in range(len(self._pos)):
            pos[i, :, :] = self._pos[i]

        # To optimize without wasting too much memory (we really have
        # troubles here) we group particles in blocks and tabulate the
        # exponentials over time this is more memory consuming but we
        # can optimize the inner loop.  even better we could change
        # order in the tabulated expo array to speed things up shape
        # is (Npart, Ndim)
        block = min(200, self._pos[0].shape[0])
        kmax = max(self.kvec.keys()) + self.dk
        acf = [defaultdict(float) for k in self.k_sorted]
        cnt = [defaultdict(float) for k in self.k_sorted]
        if self.trajectory.block_size > 1:
            skip = self.trajectory.block_size
        else:
            skip = self.skip

        for j in range(0, pos.shape[1], block):
            x = expo_sphere(self.k0, kmax, pos[:, j:j+block, :])
            for kk, knorm in enumerate(self.k_sorted):
                # Pick up a random, unique set of nk vectors out ot the avilable ones
                # without exceeding maximum number of vectors in shell nkmax
                # TODO: refactor this using _k_decimate()
                nkmax = len(self.kvec[knorm])
                for kkk in self.k_selected[kk]:
                    ik = self.kvec[knorm][kkk]
                    for off, i in self._discrete_tgrid:
                        for i0 in xrange(off, x.shape[0]-i, skip):
                            # Get the actual time difference. steps must be accessed efficiently (cached!)
                            # TODO: fix x.shape[0] in loop and x.shape[1] in normalization everywhere!
                            dt = self.trajectory.steps[i0+i] - self.trajectory.steps[i0]
                            acf[kk][dt] += numpy.sum(x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate() *
                                                     x[i0+i, :, 1, ik[1]]*x[i0, :, 1, ik[1]].conjugate() *
                                                     x[i0+i, :, 2, ik[2]]*x[i0, :, 2, ik[2]].conjugate()).real
                            cnt[kk][dt] += x.shape[1]

        t_sorted = sorted(acf[0].keys())
        self.grid[0] = self.k_sorted
        self.grid[1] = [ti*self.trajectory.timestep for ti in t_sorted]
        self.value = [[acf[kk][ti] / cnt[kk][ti] for ti in t_sorted] for kk in range(len(self.grid[0]))]
        self.value = [[self.value[kk][i] / self.value[kk][0] for i in range(len(self.value[kk]))] for kk in range(len(self.grid[0]))]

    def analyze(self):
        try:
            from .helpers import feqc
        except ImportError:
            return

        self.tau = {}
        for i, k in enumerate(self.grid[0]):
            try:
                self.tau[k] = feqc(self.grid[1], self.value[i], 1/numpy.exp(1.0))[0]
            except ValueError:
                self.tau[k] = None

    def write(self):
        Correlation.write(self)
        # TODO: refactor
        filename = '.'.join([e for e in [self.trajectory.filename, 'pp', self.short_name, self.tag] if len(e)>0])
        fileinfo = filename + '.tau'
        if not self.output is sys.stdout:
            out = open(fileinfo, 'w')
        else:
            out = sys.stdout

        # some header
        # custom writing of taus (could be refactored)
        for k in self.tau:
            if self.tau[k] is None:
                out.write('%12g\n' % k)
            else:
                out.write('%12g %12g\n' % (k, self.tau[k]))

        if not self.output is sys.stdout:
            out.close()


class IntermediateScattering(FourierSpaceCorrelation):

    nbodies = 2

    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=100, dk=0.1, tsamples=60,
                 kmin=1.0, kmax=10.0, ksamples=10):
        FourierSpaceCorrelation.__init__(self, trajectory, [kgrid, tgrid], ('k', 't'),
                                         'fkt.total', 'Intermediate scattering function',
                                         'pos', nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        check_block_size(self.trajectory.steps, self.trajectory.block_size)
        if tgrid is None:
            self.grid[1] = logx_grid(0.0, trajectory.total_time * 0.75, tsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid[1])

    def _tabulate_rho(self, k_sorted, k_selected, f=numpy.sum):

        """Tabulate densities"""

        nsteps = len(self._pos_0)
        kmax = max(self.kvec.keys()) + self.dk
        rho_0 = [defaultdict(complex) for it in range(nsteps)]
        rho_1 = [defaultdict(complex) for it in range(nsteps)]
        for it in range(nsteps):
            expo_0 = expo_sphere(self.k0, kmax, self._pos_0[it])
            # Optimize a bit here: if there is only one filter (alpha-alpha or total calculation)
            # expo_2 will be just a reference to expo_1
            if self._pos_1 is self._pos_0:
                expo_1 = expo_0
            else:
                expo_1 = expo_sphere(self.k0, kmax, self._pos_1[it])

            # Tabulate densities rho_0, rho_1
            for kk, knorm in enumerate(k_sorted):
                for i in k_selected[kk]:
                    ik = self.kvec[knorm][i]
                    rho_0[it][ik] = numpy.sum(expo_0[..., 0, ik[0]] * expo_0[..., 1, ik[1]] * expo_0[..., 2, ik[2]])
                    # Same optimization as above: only calculate rho_1 if needed
                    if not self._pos_1 is self._pos_0:
                        rho_1[it][ik] = numpy.sum(expo_1[..., 0, ik[0]] * expo_1[..., 1, ik[1]] * expo_1[..., 2, ik[2]])
            # Optimization
            if self._pos_1 is self._pos_0:
                rho_1 = rho_0

        return rho_0, rho_1

    def _compute(self):
        # Setup k vectors and tabulate densities
        k_sorted, k_selected = self._decimate_k()
        rho_0, rho_1 = self._tabulate_rho(k_sorted, k_selected)

        # Compute correlation function
        acf = [defaultdict(float) for k in k_sorted]
        cnt = [defaultdict(float) for k in k_sorted]
        skip = self.trajectory.block_size
        for kk, knorm in enumerate(k_sorted):
            for j in k_selected[kk]:
                ik = self.kvec[knorm][j]
                for off, i in self._discrete_tgrid:
                    for i0 in xrange(off, len(rho_0)-i, skip):
                        # Get the actual time difference
                        # TODO: It looks like the order of i0 and ik lopps should be swapped
                        dt = self.trajectory.steps[i0+i] - self.trajectory.steps[i0]
                        acf[kk][dt] += (rho_0[i0+i][ik] * rho_1[i0][ik].conjugate()).real #/ self._pos[i0].shape[0]
                        cnt[kk][dt] += 1

        # Normalization
        times = sorted(acf[0].keys())
        self.grid[0] = k_sorted
        self.grid[1] = [ti*self.trajectory.timestep for ti in times]
        if self._pos_0 is self._pos_1:
            # First normalize by cnt (time counts), then by value at t=0
            # We do not need to normalize by the average number of particles
            # TODO: check normalization when not GC, does not give exactly the short time behavior as pp.x
            nav = sum([p.shape[0] for p in self._pos]) / len(self._pos)
            self.value_nonorm = [[acf[kk][ti] / (cnt[kk][ti]) for ti in times] for kk in range(len(self.grid[0]))]
            self.value = [[v / self.value_nonorm[kk][0] for v in self.value_nonorm[kk]] for kk in range(len(self.grid[0]))]
        else:
            nav_0 = sum([p.shape[0] for p in self._pos_0]) / len(self._pos_0)
            nav_1 = sum([p.shape[0] for p in self._pos_1]) / len(self._pos_1)
            self.value_nonorm = [[acf[kk][ti] / (cnt[kk][ti]) for ti in times] for kk in range(len(self.grid[0]))]
            self.value = [[v / self.value_nonorm[kk][0] for v in self.value_nonorm[kk]] for kk in range(len(self.grid[0]))]

    def analyze(self):
        try:
            from .helpers import feqc
        except ImportError:
            return
        self.tau = {}
        for i, k in enumerate(self.grid[0]):
            try:
                self.tau[k] = feqc(self.grid[1], self.value[i], 1/numpy.exp(1.0))[0]
            except ValueError:
                self.tau[k] = None

    def write(self):
        Correlation.write(self)
        # Write down unnormalized functions
        # Correlation.write(self, self.value_nonorm)

        # TODO: refactor
        filename = '.'.join([e for e in [self.trajectory.filename, 'pp', self.short_name, self.tag] if len(e)>0])
        fileinfo = filename + '.tau'
        if not self.output is sys.stdout:
            out = open(fileinfo, 'w')
        else:
            out = sys.stdout

        # some header
        # custom writing of taus (could be refactored)
        for k in self.tau:
            if self.tau[k] is None:
                out.write('%12g\n' % k)
            else:
                out.write('%12g %12g\n' % (k, self.tau[k]))

        if not self.output is sys.stdout:
            out.close()


class StructureFactor(FourierSpaceCorrelation):

    nbodies = 2

    def __init__(self, trajectory, kgrid=None, norigins=-1, nk=20,
                 dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30,
                 trajectory_field=None, field=None):
        """
        If `trajectory_field` is not None, the field is read from the last
        column of this trajectory file, unless the `field` string is
        provided.
        """
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, 'k',
                                         'sk', 'structure factor',
                                         ['pos'], nk, dk, kmin,
                                         kmax, ksamples)
        # TODO: move this up the chain?
        self.skip = adjust_skip(self.trajectory, norigins)
        self._is_cell_variable = None
        self._field, tag = self._add_field(trajectory_field, field)
        self.tag += tag

    def _add_field(self, trajectory_field, field):
        if trajectory_field is None:
            return None, ''
        else:
            # TODO: check step consistency 06.09.2017
            from atooms.trajectory import TrajectoryXYZ
            with TrajectoryXYZ(trajectory_field) as th:
                fields = []
                # This must be a string, not a list
                unique_field = th._read_metadata(0)['columns']
                if isinstance(unique_field, list):
                    if field is None:
                        unique_field = unique_field[-1]
                    else:
                        unique_field = field
                for s in th:
                    fields.append(s.dump('particle.%s' % unique_field))
            return fields, unique_field

    def _compute(self):
        from atooms.trajectory.utils import is_cell_variable
        nsteps = len(self._pos_0)
        # Setup k vectors and tabulate rho
        k_sorted, k_selected = self._decimate_k()
        kmax = max(self.kvec.keys()) + self.dk
        cnt = [0 for k in k_sorted]
        rho_av = [complex(0.,0.) for k in k_sorted]
        rho2_av = [complex(0.,0.) for k in k_sorted]
        variable_cell = is_cell_variable(self.trajectory)
        for i in range(0, nsteps, self.skip):
            # If cell changes we have to update the wave vectors
            if variable_cell:
                self._setup(i)
                k_sorted, k_selected = self._decimate_k()
                kmax = max(self.kvec.keys()) + self.dk

            # Tabulate exponentials
            if self._pos_0[i] is self._pos_1[i]:
                # Identical species
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_0
            else:
                # Cross correlation
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_sphere(self.k0, kmax, self._pos_1[i])

            for kk, knorm in enumerate(k_sorted):
                for k in k_selected[kk]:
                    ik = self.kvec[knorm][k]
                    # In the absence of a microscopic field, rho_av = (0, 0)
                    if not self._field:
                        if expo_0 is expo_1:
                            # Identical species
                            rho_0 = numpy.sum(expo_0[...,0,ik[0]] *
                                              expo_0[...,1,ik[1]] *
                                              expo_0[...,2,ik[2]])
                            rho_1 = rho_0
                        else:
                            # Cross correlation
                            rho_0 = numpy.sum(expo_0[...,0,ik[0]] *
                                              expo_0[...,1,ik[1]] *
                                              expo_0[...,2,ik[2]])
                            rho_1 = numpy.sum(expo_1[...,0,ik[0]] *
                                              expo_1[...,1,ik[1]] *
                                              expo_1[...,2,ik[2]])
                    else:
                        # We have a field as a weight
                        rho_0 = numpy.sum(self._field[i] *
                                          expo_0[...,0,ik[0]] *
                                          expo_0[...,1,ik[1]] *
                                          expo_0[...,2,ik[2]])
                        rho_1 = rho_0
                        rho_av[kk] += rho_0

                    rho2_av[kk] += (rho_0 * rho_1.conjugate())
                    cnt[kk] += 1

        # Normalization.
        npart_0 = sum([p.shape[0] for p in self._pos_0]) / float(len(self._pos_0))
        npart_1 = sum([p.shape[0] for p in self._pos_1]) / float(len(self._pos_1))
        self.grid = k_sorted
        self.value, self.value_nonorm = [], []
        for kk in range(len(self.grid)):
            norm = float(npart_0 * npart_1)**0.5
            value = (rho2_av[kk] / cnt[kk] - rho_av[kk]*rho_av[kk].conjugate() / cnt[kk]**2).real            
            self.value.append(value / norm)
            self.value_nonorm.append(value)


class SpectralDensity(FourierSpaceCorrelation):

    """
    Free volume spectral density.

    See Zachary, Jiao, Torquato PRL 106, 178001 (2011).
    """

    def __init__(self, trajectory, trajectory_radius, kgrid=None,
                 norigins=-1, nk=20, dk=0.1, kmin=-1.0, kmax=15.0,
                 ksamples=30):
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, 'k',
                                         'ik', 'spectral density',
                                         ['pos'], nk, dk, kmin,
                                         kmax, ksamples)
        # TODO: move this up the chain?
        self.skip = adjust_skip(self.trajectory, norigins)
        self._is_cell_variable = None
        # TODO: check step consistency 06.09.2017
        from atooms.trajectory import TrajectoryXYZ, Trajectory
        with Trajectory(trajectory_radius) as th:
            self._radius = [s.dump('particle.radius') for s in th]

    def _compute(self):
        nsteps = len(self._pos)
        # Setup k vectors and tabulate rho
        k_sorted, k_selected = self._decimate_k()
        kmax = max(self.kvec.keys()) + self.dk
        cnt = [0 for k in k_sorted]
        # Note: actually rho_av is not calculated because it is negligible
        rho_av = [complex(0.,0.) for k in k_sorted]
        rho2_av = [complex(0.,0.) for k in k_sorted]
        cell_variable = is_cell_variable(self.trajectory)
        for i in range(0, nsteps, self.skip):
            # If cell changes we have to update
            if cell_variable:
                self._setup(i)
                k_sorted, k_selected = self._decimate_k()
                kmax = max(self.kvec.keys()) + self.dk

            expo = expo_sphere(self.k0, kmax, self._pos[i])
            for kk, knorm in enumerate(k_sorted):
                for k in k_selected[kk]:
                    ik = self.kvec[knorm][k]
                    Ri = self._radius[i]
                    mk = 4 * numpy.pi / knorm**3 * (numpy.sin(knorm*Ri) - (knorm*Ri) * numpy.cos(knorm*Ri))
                    rho = numpy.sum(mk*expo[...,0,ik[0]]*expo[...,1,ik[1]]*expo[...,2,ik[2]])
                    rho2_av[kk] += (rho * rho.conjugate())
                    cnt[kk] += 1

        # Normalization.
        volume = numpy.average([s.cell.volume for s in self.trajectory])
        self.grid = k_sorted
        self.value = [(rho2_av[kk] / cnt[kk] - rho_av[kk]*rho_av[kk].conjugate() / cnt[kk]**2).real / volume
                       for kk in range(len(self.grid))]
        self.value_nonorm = [rho2_av[kk].real / cnt[kk]
                             for kk in range(len(self.grid))]


class StructureFactorStats(FourierSpaceCorrelation):

    def __init__(self, trajectory, kgrid=None, norigins=-1, nk=1000, dk=1.0, kmin=7.0):
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, 'k', 'skstats', 'structure factor stats', ['pos'], \
                                         nk, dk, kmin, kmin, 1)
        # TODO: move this up the chain?
        self.skip = adjust_skip(self.trajectory, norigins)

    def _compute(self):
        def skew(x):
            return numpy.sum((x-numpy.mean(x))**3) / len(x) / numpy.std(x)**3

        # Setup k vectors and tabulate rho
        k_sorted, k_selected = self._decimate_k()
        nsteps = len(self._pos)
        kmax = max(self.kvec.keys()) + self.dk
        self._mean = []; self._var = []; self._skew = []
        self.grid = range(0, nsteps, self.skip)
        for i in range(0, nsteps, self.skip):
            cnt = 0
            sk = []
            expo = expo_sphere(self.k0, kmax, self._pos[i])
            npart = self._pos[i].shape[0]
            for kk, knorm in enumerate(k_sorted):
                for k in k_selected[kk]:
                    ik = self.kvec[knorm][k]
                    rho = numpy.sum(expo[...,0,ik[0]] * expo[...,1,ik[1]] * expo[...,2,ik[2]])
                    rho2 = rho * rho.conjugate()
                    sk.append(rho2.real / npart)
            self._mean.append(numpy.average(sk))
            self._var.append(numpy.var(sk))
            self._skew.append(skew(sk))

    def write(self):
        comments = """\
# ave = %g
# var = %g
# skew = %g
""" % (numpy.average(self._mean), numpy.average(self._var), numpy.average(self._skew))
        with open(self._output_file, 'w') as fh:
            fh.write(comments)
            numpy.savetxt(fh, numpy.array(zip(self.grid, self._var, self._skew)), fmt="%g")


class S4ktOverlap(FourierSpaceCorrelation):

    # TODO: refactor a S4k base correlation that forces to implement tabulat method (e.g. overlap, Q_6, voronoi ...)
    # TODO: should we drop this instead and rely on F(k,t) with grandcanonical

    def __init__(self, trajectory, tgrid, kgrid=None, norigins=-1, nk=20, dk=0.1, a=0.3, kmin=1.0, kmax=10.0, ksamples=10):
        FourierSpaceCorrelation.__init__(self, trajectory, [tgrid, kgrid], ('t', 'k'), 's4kt', \
                                         '4-point dynamic structure factor', ['pos','pos-unf'], \
                                         nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        self._discrete_tgrid = setup_t_grid(trajectory, tgrid)
        self.skip = adjust_skip(self.trajectory, norigins)
        self.a_square = a**2

        # Setup k vectors and tabulate rho
        # TODO: move decimate up the chain?
        self.k_sorted, self.k_selected = self._decimate_k()
        # Redefine kgrid to give exactly the average wave vectors used.
        # TODO; should we do it for in base?
        self.grid[1] = self._actual_k_grid(self.k_sorted, self.k_selected)

    def _tabulate_W(self, k_sorted, k_selected, t_off, t, skip):
        """ Tabulate W """
        nsteps = len(self._pos)
        side = self.trajectory[0].cell.side
        kmax = max(self.kvec.keys()) + self.dk
        nt = xrange(t_off, len(self._pos)-t, skip)
        W = {}
        for i_0, t_0 in enumerate(nt):
            expo = expo_sphere(self.k0, kmax, self._pos[t_0])
            for kk, knorm in enumerate(k_sorted):
                for i in k_selected[kk]:
                    ik = self.kvec[knorm][i]
                    if not ik in W:
                        W[ik] = numpy.ndarray(len(nt), dtype=complex)
                    W[ik][i_0] = numpy.sum(self_overlap(self._pos_unf[t_0], self._pos_unf[t_0+t], side, self.a_square) *
                                          expo[...,0,ik[0]] * expo[...,1,ik[1]] * expo[...,2,ik[2]])
        return W

    def _compute(self):
        # Make sure there is only one time in tgrid.
        # We could easily workaround it by outer looping over i
        # We do not expect to do it for many times (typically we show S_4(k,tau_alpha) vs k)
        # if len(self._discrete_tgrid) > 1:
        #     raise ValueError('There should be only one time for S4kt')
        dt = []
        self.value = []
        for off, i  in self._discrete_tgrid:

            # as for fkt
            W = self._tabulate_W(self.k_sorted, self.k_selected, off, i, self.skip)

            # Compute vriance of W
            cnt = [0 for k in self.k_sorted]
            w_av = [complex(0., 0.) for k in self.k_sorted]
            w2_av = [complex(0., 0.) for k in self.k_sorted]
            for kk, knorm in enumerate(self.k_sorted):
                for j in self.k_selected[kk]:
                    ik = self.kvec[knorm][j]
                    # Comupte |<W>|^2  and <W W^*>
                    w_av[kk] = numpy.average(W[ik])
                    w2_av[kk] = numpy.average(W[ik] * W[ik].conjugate())

            # Normalization
            npart = self._pos[0].shape[0]
            dt.append(self.trajectory.timestep * (self.trajectory.steps[off+i] - self.trajectory.steps[off]))
            #self.grid[1] = k_sorted
            self.value.append([float(w2_av[kk] - (w_av[kk]*w_av[kk].conjugate())) / npart for kk in range(len(self.grid[1]))])
        self.grid[0] = dt
