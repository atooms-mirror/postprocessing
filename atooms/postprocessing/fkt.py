# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Intermediate scattering function."""

from collections import defaultdict

import numpy
from atooms.trajectory.utils import check_block_size

from .helpers import logx_grid, setup_t_grid
from .correlation import Correlation
from .fourierspace import FourierSpaceCorrelation, expo_sphere
from .progress import progress

__all__ = ['SelfIntermediateScattering', 'IntermediateScattering']


def _write_tau(out, db):
    # Custom writing of relaxation times
    out.write('# title: relaxation times tau(k) as a function of k\n')
    out.write('# columns: k, tau(k)\n')
    out.write('# note: tau is the time at which the correlation function has decayed to 1/e\n')
    for k, tau in db['relaxation times tau'].items():
        if tau is None:
            out.write('%g\n' % k)
        else:
            out.write('%g %g\n' % (k, tau))

def _extract_tau(k, t, f):
    from .helpers import feqc
    tau = {}
    for i, k in enumerate(k):
        try:
            tau[k] = feqc(t, f[i], 1 / numpy.exp(1.0))[0]
        except ValueError:
            tau[k] = None
    return tau


class SelfIntermediateScattering(FourierSpaceCorrelation):
    """
    Self part of the intermediate scattering function.

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    symbol = 'fskt'
    short_name = 'F_s(k,t)'
    long_name = 'self intermediate scattering function'
    phasespace = 'pos-unf'

    #TODO: xyz files are 2 slower than hdf5 here
    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=8, tsamples=60,
                 dk=0.1, kmin=1.0, kmax=10.0, ksamples=10, norigins=-1):
        FourierSpaceCorrelation.__init__(self, trajectory, [kgrid, tgrid], norigins,
                                         nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        # Before setting up the time grid, we need to check periodicity over blocks
        check_block_size(self.trajectory.steps, self.trajectory.block_size)
        if tgrid is None:
            self.grid[1] = [0.0] + logx_grid(self.trajectory.timestep,
                                             self.trajectory.total_time * 0.75, tsamples)
        self._discrete_tgrid = setup_t_grid(self.trajectory, self.grid[1])

    def _compute(self):
        # Throw everything into a big numpy array (nframes, npos, ndim)
        pos = numpy.array(self._pos_unf)

        # To optimize without wasting too much memory (we have
        # troubles here) we group particles in blocks and tabulate the
        # exponentials over time. This is more memory consuming but we
        # can optimize the inner loop. Even better, we could change
        # the order in the tabulated expo array to speed things up
        # Use 10 blocks, but do not exceed 200 particles
        number_of_blocks = 10
        block = int(self._pos_unf[0].shape[0] / float(number_of_blocks))
        block = max(20, block)
        block = min(200, block)
        if len(self.kvector.keys()) == 0:
            raise ValueError('could not find any wave-vectors, try increasing dk')
        kmax = max(self.kvector.keys()) + self.dk
        acf = [defaultdict(float) for _ in self.kgrid]
        cnt = [defaultdict(float) for _ in self.kgrid]
        if self.trajectory.block_size > 1:
            skip = self.trajectory.block_size
        else:
            skip = self.skip

        origins = range(0, pos.shape[1], block)
        for j in progress(origins):
            x = expo_sphere(self.k0, kmax, pos[:, j:j + block, :])
            for kk, knorm in enumerate(self.kgrid):
                for kkk in self.selection[kk]:
                    ik = self.kvector[knorm][kkk]
                    for off, i in self._discrete_tgrid:
                        for i0 in range(off, x.shape[0]-i, skip):
                            # Get the actual time difference. steps must be accessed efficiently (cached!)
                            dt = self.trajectory.steps[i0+i] - self.trajectory.steps[i0]
                            acf[kk][dt] += numpy.sum(x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate() *
                                                     x[i0+i, :, 1, ik[1]]*x[i0, :, 1, ik[1]].conjugate() *
                                                     x[i0+i, :, 2, ik[2]]*x[i0, :, 2, ik[2]].conjugate()).real
                            cnt[kk][dt] += x.shape[1]

        tgrid = sorted(acf[0].keys())
        self.grid[0] = self.kgrid
        self.grid[1] = [ti*self.trajectory.timestep for ti in tgrid]
        self.value = [[acf[kk][ti] / cnt[kk][ti] for ti in tgrid] for kk in range(len(self.grid[0]))]
        self.value = [[self.value[kk][i] / self.value[kk][0] for i in range(len(self.value[kk]))] for kk in range(len(self.grid[0]))]

    def analyze(self):
        self.analysis['relaxation times tau'] = _extract_tau(self.grid[0], self.grid[1], self.value)

    def write(self):
        Correlation.write(self)
        if self._output_file != '/dev/stdout':
            with open(self._output_file + '.tau', 'w') as out:
                _write_tau(out, self.analysis)


class IntermediateScattering(FourierSpaceCorrelation):
    """
    Coherent intermediate scattering function.

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    nbodies = 2
    symbol = 'fkt'
    short_name = 'F(k,t)'
    long_name = 'intermediate scattering function'
    phasespace = 'pos'

    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=100, dk=0.1, tsamples=60,
                 kmin=1.0, kmax=10.0, ksamples=10, norigins=-1):
        FourierSpaceCorrelation.__init__(self, trajectory, [kgrid, tgrid], norigins,
                                         nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        check_block_size(self.trajectory.steps, self.trajectory.block_size)
        if tgrid is None:
            self.grid[1] = logx_grid(0.0, self.trajectory.total_time * 0.75, tsamples)
        self._discrete_tgrid = setup_t_grid(self.trajectory, self.grid[1])

    def _tabulate_rho(self, kgrid, selection):
        """
        Tabulate densities
        """
        nsteps = len(self._pos_0)
        if len(self.kvector.keys()) == 0:
            raise ValueError('could not find any wave-vectors, try increasing dk')
        kmax = max(self.kvector.keys()) + self.dk
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
            for kk, knorm in enumerate(kgrid):
                for i in selection[kk]:
                    ik = self.kvector[knorm][i]
                    rho_0[it][ik] = numpy.sum(expo_0[..., 0, ik[0]] * expo_0[..., 1, ik[1]] * expo_0[..., 2, ik[2]])
                    # Same optimization as above: only calculate rho_1 if needed
                    if self._pos_1 is not self._pos_0:
                        rho_1[it][ik] = numpy.sum(expo_1[..., 0, ik[0]] * expo_1[..., 1, ik[1]] * expo_1[..., 2, ik[2]])
            # Optimization
            if self._pos_1 is self._pos_0:
                rho_1 = rho_0

        return rho_0, rho_1

    def _compute(self):
        # Setup k vectors and tabulate densities
        kgrid, selection =  self.kgrid, self.selection
        rho_0, rho_1 = self._tabulate_rho(kgrid, selection)

        # Compute correlation function
        acf = [defaultdict(float) for _ in kgrid]
        cnt = [defaultdict(float) for _ in kgrid]
        skip = self.trajectory.block_size
        for kk, knorm in enumerate(progress(kgrid)):
            for j in selection[kk]:
                ik = self.kvector[knorm][j]
                for off, i in self._discrete_tgrid:
                    for i0 in range(off, len(rho_0)-i, skip):
                        # Get the actual time difference
                        # TODO: It looks like the order of i0 and ik lopps should be swapped
                        dt = self.trajectory.steps[i0+i] - self.trajectory.steps[i0]
                        acf[kk][dt] += (rho_0[i0+i][ik] * rho_1[i0][ik].conjugate()).real #/ self._pos[i0].shape[0]
                        cnt[kk][dt] += 1

        # Normalization
        times = sorted(acf[0].keys())
        self.grid[0] = kgrid
        self.grid[1] = [ti*self.trajectory.timestep for ti in times]
        if self._pos_0 is self._pos_1:
            # First normalize by cnt (time counts), then by value at t=0
            # We do not need to normalize by the average number of particles
            # TODO: check normalization when not GC, does not give exactly the short time behavior as pp.x
            nav = sum([p.shape[0] for p in self._pos]) / len(self._pos)
            self.value_nonorm = [[acf[kk][ti] / (cnt[kk][ti]) for ti in times] for kk in range(len(self.grid[0]))]
            self.value = [[v / self.value_nonorm[kk][0] for v in self.value_nonorm[kk]] for kk in range(len(self.grid[0]))]
        else:
            # nav_0 = sum([p.shape[0] for p in self._pos_0]) / len(self._pos_0)
            # nav_1 = sum([p.shape[0] for p in self._pos_1]) / len(self._pos_1)
            self.value_nonorm = [[acf[kk][ti] / (cnt[kk][ti]) for ti in times] for kk in range(len(self.grid[0]))]
            self.value = [[v / self.value_nonorm[kk][0] for v in self.value_nonorm[kk]] for kk in range(len(self.grid[0]))]

    def analyze(self):
        self.analysis['relaxation times tau'] = _extract_tau(self.grid[0], self.grid[1], self.value)

    def write(self):
        Correlation.write(self)
        if self._output_file != '/dev/stdout':
            with open(self._output_file + '.tau', 'w') as out:
                _write_tau(out, self.analysis)
