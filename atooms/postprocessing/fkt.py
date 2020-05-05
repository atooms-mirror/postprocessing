# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Intermediate scattering function."""

import logging
from collections import defaultdict

import numpy
from atooms.trajectory.utils import check_block_size

from . import core
from atooms.trajectory import Trajectory
from .helpers import logx_grid, setup_t_grid
from .correlation import Correlation
from .fourierspace import FourierSpaceCorrelation, expo_sphere
from .progress import progress

__all__ = ['SelfIntermediateScattering',
           'SelfIntermediateScatteringLegacy', 'SelfIntermediateScatteringFast',
           'IntermediateScattering']

_log = logging.getLogger(__name__)


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

def _extract_tau(k, t, f, factor=1/numpy.exp(1.0)):
    from .helpers import feqc
    # Ensure first point in time grid is t=0
    # This will be enforced by the classes below but we never know
    if t[0] > 0:
        raise ValueError('First point in time grid must be zero')
    tau = {}
    for i, k in enumerate(k):
        try:
            tau[k] = feqc(t, f[i], f[i][0]*factor)[0]
        except ValueError:
            tau[k] = None
    return tau


class IntermediateScatteringBase(FourierSpaceCorrelation):

    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=1,
                 tsamples=1, dk=0.1, kmin=1.0, kmax=10.0, ksamples=10,
                 norigins=-1, fix_cm=False, normalize=True):
        FourierSpaceCorrelation.__init__(self, trajectory, [kgrid, tgrid], norigins,
                                         nk, dk, kmin, kmax, ksamples, fix_cm, normalize)
        # Before setting up the time grid, we need to check periodicity over blocks
        try:
            check_block_size(self.trajectory.steps, self.trajectory.block_size)
        except IndexError as e:
            _log.warn('issue with trajectory blocks, the time grid may not correspond to the requested one (%s)', e.message)

        # Setup time grid
        if self.grid[1] is None:
            self.grid[1] = logx_grid(0.0, self.trajectory.total_time * 0.75, tsamples)
        else:
            # If the values are normalized, we make sure the
            # user-provided time grid includes t=0. It is removed
            # after normalization
            if self.normalize:
                if self.grid[1][0] > 0:
                    _log.info('adding t=0 to the time grid to normalize F_s(k,t)')
                    self.grid[1] = [0.0] + self.grid[1]

        # When a single time origin is requested,
        # make sure no additional time origins except the first frame is used
        self._discrete_tgrid = setup_t_grid(self.trajectory, self.grid[1], offset=norigins != '1' and norigins != 1)


class SelfIntermediateScatteringLegacy(IntermediateScatteringBase):
    """
    Self part of the intermediate scattering function.

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    symbol = 'fskt'
    short_name = 'F_s(k,t)'
    long_name = 'self intermediate scattering function'
    phasespace = 'pos'

    def __init__(self, trajectory, kgrid=None, tgrid=None, nk=8,
                 tsamples=60, dk=0.1, kmin=1.0, kmax=10.0,
                 ksamples=10, norigins=-1, fix_cm=False,
                 lookup_mb=64.0, normalize=True):
        # TODO: remove this backward compatible tgrid fix in a major release
        # The default time grid should be the same for F_s(k,t) and F(k,t)
        if isinstance(trajectory, str):
            trajectory = Trajectory(trajectory, mode='r', fmt=core.pp_trajectory_format)
        if tgrid is None:
            tgrid = [0.0] + logx_grid(trajectory.timestep,
                                      trajectory.total_time * 0.75, tsamples)
        super(SelfIntermediateScatteringLegacy,
              self).__init__(trajectory, kgrid=kgrid, tgrid=tgrid,
                             nk=nk, tsamples=tsamples, dk=dk, kmin=kmin,
                             kmax=kmax, ksamples=ksamples, norigins=norigins,
                             fix_cm=fix_cm, normalize=normalize)
        self.lookup_mb = lookup_mb
        """Memory in Mb allocated for exponentials tabulation"""

    def _compute(self):
        # Throw everything into a big numpy array (nframes, npos, ndim)
        pos = numpy.array(self._pos)
        ndims = len(self.k0)
        # To optimize without wasting too much memory (we have
        # troubles here) we group particles in blocks and tabulate the
        # exponentials over time. This is more memory consuming but we
        # can optimize the inner loop. Even better, we could change
        # the order in the tabulated expo array to speed things up
        # Use 10 blocks, but do not exceed 200 particles
        number_of_blocks = 10
        block = int(pos[0].shape[0] / float(number_of_blocks))
        block = max(20, block)
        block = min(200, block)
        if len(self.kvector.keys()) == 0:
            raise ValueError('could not find any wave-vectors, try increasing dk')
        kmax = max(self.kvector.keys()) + self.dk
        acf = [defaultdict(float) for _ in self.kgrid]
        cnt = [defaultdict(float) for _ in self.kgrid]
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
                            # Dimensional switch
                            if ndims == 3:
                                acf[kk][dt] += numpy.sum(x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate() *
                                                         x[i0+i, :, 1, ik[1]]*x[i0, :, 1, ik[1]].conjugate() *
                                                         x[i0+i, :, 2, ik[2]]*x[i0, :, 2, ik[2]].conjugate()).real
                            elif ndims == 2:
                                acf[kk][dt] += numpy.sum(x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate() *
                                                         x[i0+i, :, 1, ik[1]]*x[i0, :, 1, ik[1]].conjugate()).real

                            else:
                                # Arbitrary dimension (a bit slower)
                                tmp = x[i0+i, :, 0, ik[0]]*x[i0, :, 0, ik[0]].conjugate()
                                for idim in range(1, len(ik)):
                                    tmp *= x[i0+i, :, idim, ik[idim]]*x[i0, :, idim, ik[idim]].conjugate()
                                acf[kk][dt] += numpy.sum(tmp).real                            
                            cnt[kk][dt] += x.shape[1]

        tgrid = sorted(acf[0].keys())
        self.grid[0] = self.kgrid
        self.grid[1] = [ti*self.trajectory.timestep for ti in tgrid]
        self.value = [[acf[kk][ti] / cnt[kk][ti] for ti in tgrid] for kk in range(len(self.grid[0]))]

        # Normalize
        if self.normalize:
            for k in range(len(self.grid[0])):
                for i in range(len(self.value[k])):
                    self.value[k][i] /= self.value[k][0]

    def analyze(self):
        self.analysis['relaxation times tau'] = _extract_tau(self.grid[0], self.grid[1], self.value)

    def write(self):
        Correlation.write(self)
        if self._output_file != '/dev/stdout':
            with open(self._output_file + '.tau', 'w') as out:
                _write_tau(out, self.analysis)


class SelfIntermediateScatteringFast(SelfIntermediateScatteringLegacy):
    """
    Self part of the intermediate scattering function (fast version)
    
    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """        
    def _compute(self):
        try:
            from atooms.postprocessing.fourierspace_wrap import fourierspace_module
        except ImportError:
            _log.error('f90 wrapper missing or not functioning')
            raise

        # Throw everything into a big numpy array (nframes, npos, ndim)
        pos = numpy.array(self._pos)

        # Select the f90 kernel
        ndims = len(self.k0)
        if ndims == 3:
            fskt_kernel = fourierspace_module.fskt_kernel_3d
        elif ndims == 2:
            fskt_kernel = fourierspace_module.fskt_kernel_2d
        else:
            fskt_kernel = fourierspace_module.fskt_kernel_nd
            
        # To optimize without wasting too much memory (we have
        # troubles here) we group particles in blocks and tabulate the
        # exponentials over time. This is more memory consuming but we
        # can optimize the inner loop. The esitmated amuount of
        # allocated memory in Mb for the expo array is
        # self.lookup_mb. Note that the actual memory need scales
        # with number of k vectors, system size and number of frames.
        kmax = max(self.kvector.keys()) + self.dk
        kvec_size = 2*(1 + int(kmax / min(self.k0))) + 1
        pos_size = numpy.product(pos.shape)
        target_size = self.lookup_mb * 1e6 / 16.  # 16 bytes for a (double) complex        
        number_of_blocks = int(pos_size * kvec_size / target_size)
        number_of_blocks = max(1, number_of_blocks)
        block = int(pos[0].shape[0] / float(number_of_blocks))
        block = max(1, block)
        block = min(pos.shape[1], block)
        if len(self.kvector.keys()) == 0:
            raise ValueError('could not find any wave-vectors, try increasing dk')
        acf = [defaultdict(float) for _ in self.kgrid]
        cnt = [defaultdict(float) for _ in self.kgrid]
        skip = self.skip
        origins = range(0, pos.shape[1], block)
        for j in progress(origins):
            x = expo_sphere(self.k0, kmax, pos[:, j:j + block, :])            
            xf = numpy.asfortranarray(x)
            for kk, knorm in enumerate(self.kgrid):
                for kkk in self.selection[kk]:
                    ik = self.kvector[knorm][kkk]
                    for off, i in self._discrete_tgrid:
                        for i0 in range(off, x.shape[0]-i, skip):
                            # Get the actual time difference. steps must be accessed efficiently (cached!)
                            dt = self.trajectory.steps[i0+i] - self.trajectory.steps[i0]
                            # Call f90 kernel
                            res = fskt_kernel(xf, i0+1, i0+1+i, numpy.array(ik, dtype=numpy.int32)+1)
                            acf[kk][dt] += res.real
                            cnt[kk][dt] += x.shape[1]
                            
        tgrid = sorted(acf[0].keys())
        self.grid[0] = self.kgrid
        self.grid[1] = [ti*self.trajectory.timestep for ti in tgrid]
        self.value = [[acf[kk][ti] / cnt[kk][ti] for ti in tgrid] for kk in range(len(self.grid[0]))]
        # Normalize
        if self.normalize:
            for k in range(len(self.grid[0])):
                for i in range(len(self.value[k])):
                    self.value[k][i] /= self.value[k][0]


# Defaults to fast
SelfIntermediateScattering = SelfIntermediateScatteringFast

class IntermediateScattering(IntermediateScatteringBase):
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
                 kmin=1.0, kmax=10.0, ksamples=10, norigins=-1, fix_cm=False, normalize=True):
        super(IntermediateScattering, self).__init__(trajectory, kgrid=kgrid, tgrid=tgrid,
                             nk=nk, tsamples=tsamples, dk=dk, kmin=kmin,
                             kmax=kmax, ksamples=ksamples, norigins=norigins,
                             fix_cm=fix_cm, normalize=normalize)

    def _tabulate_rho(self, kgrid, selection):
        """
        Tabulate densities
        """
        nsteps = len(self._pos_0)
        ndims = len(self.k0)
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
                    if ndims == 3:
                        rho_0[it][ik] = numpy.sum(expo_0[..., 0, ik[0]] * expo_0[..., 1, ik[1]] * expo_0[..., 2, ik[2]])
                    elif ndims == 2:
                        rho_0[it][ik] = numpy.sum(expo_0[..., 0, ik[0]] * expo_0[..., 1, ik[1]])
                    else:
                        # Arbitrary dimension (a bit slower)
                        tmp = expo_0[..., 0, ik[0]]
                        for idim in range(1, len(ik)):
                            tmp *= expo_0[..., idim, ik[idim]]
                        rho_0[it][ik] += numpy.sum(tmp).real                            
                        
                    # Same optimization as above: only calculate rho_1 if needed
                    if self._pos_1 is not self._pos_0:
                        if ndims == 3:
                            rho_1[it][ik] = numpy.sum(expo_1[..., 0, ik[0]] * expo_1[..., 1, ik[1]] * expo_1[..., 2, ik[2]])
                        elif ndims == 2:
                            rho_1[it][ik] = numpy.sum(expo_1[..., 0, ik[0]] * expo_1[..., 1, ik[1]])
                        else:
                            # Arbitrary dimension (a bit slower)
                            tmp = expo_1[..., 0, ik[0]]
                            for idim in range(1, len(ik)):
                                tmp *= expo_1[..., idim, ik[idim]]
                            rho_1[it][ik] += numpy.sum(tmp).real                            

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
        skip = self.skip
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
        # TODO: check normalization when not GC, does not give exactly the short time behavior as pp.x
        # This switch is not used
        # if self._pos_0 is self._pos_1:
        #     nav = sum([p.shape[0] for p in self._pos]) / len(self._pos)
        # else:
        #     nav_0 = sum([p.shape[0] for p in self._pos_0]) / len(self._pos_0)
        #     nav_1 = sum([p.shape[0] for p in self._pos_1]) / len(self._pos_1)
        # First normalize by cnt (time counts), then by value at t=0
        # We do not need to normalize by the average number of particles
        self.value_nonorm = [[acf[kk][ti] / (cnt[kk][ti]) for ti in times] for kk in range(len(self.grid[0]))]

        # Normalize
        if self.normalize:
            self.value = [[v / self.value_nonorm[kk][0] for v in self.value_nonorm[kk]] for kk in range(len(self.grid[0]))]
        else:
            self.value = self.value_nonorm

    def analyze(self):
        self.analysis['relaxation times tau'] = _extract_tau(self.grid[0], self.grid[1], self.value)

    def write(self):
        Correlation.write(self)
        if self._output_file != '/dev/stdout':
            with open(self._output_file + '.tau', 'w') as out:
                _write_tau(out, self.analysis)
