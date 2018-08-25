# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Mean square displacement."""

import numpy
import logging

from .helpers import linear_grid
from .correlation import Correlation, gcf_offset
from .helpers import adjust_skip, setup_t_grid

__all__ = ['MeanSquareDisplacement']

log = logging.getLogger(__name__)


def partition(inp, nbl):
    nel = len(inp) // nbl
    a = []
    for i in range(nbl):
        a.append(slice(i * nel, (i+1) * nel))
    return a


class MeanSquareDisplacement(Correlation):

    """Mean square displacement."""

    def __init__(self, trajectory, tgrid=None, sigma=1.0, norigins=50,
                 nsamples=30, sigma_max=1e100, nblocks=1):
        # TODO: optimize targeting msd takes a lot of time especially on large systems because of pbc unfolding
        # TODO: memory is leaking when sourcing and analyzing many files!
        self.sigma = sigma
        self.sigma_max = sigma_max
        self.nblocks = nblocks
        self.var = None

        Correlation.__init__(self, trajectory, tgrid, 't', 'msd',
                             'mean square displacement dr^2(t)', ['pos-unf'])
        if not self._need_update:
            return
        # TODO: subtrajectories should behave well when sampling is logarithmic

        # We redefine trajectory here to avoid unfolding the file if this is not necessary
        # e.g. because we are just sourcing the correlation object from file
        if tgrid is None:
            if sigma_max < 1e99:
                self.grid = linear_grid(0.0, min(trajectory.total_time * 1./(1+nblocks),
                                                 trajectory.time_when_msd_is(sigma_max**2)),
                                        nsamples)
            else:
                self.grid = linear_grid(0.0, trajectory.total_time * 1./(1+nblocks), nsamples)
        else:
            self.grid = tgrid
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid)
        self.skip = adjust_skip(trajectory, norigins)

    def _compute(self):
        # We could compute the individual directions (x,y,z) and
        # average them, which will save time and space and give an
        # estimate of the incertainty.
        def f_all(x, y):
            return numpy.sum((x-y)**2) / float(x.shape[0])
        def f(x, y):
            return numpy.sum((x-y)**2) / float(x.shape[0])

        # Go for it. Redefine the time grid.
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid, self.skip,
                                           self.trajectory.steps, self._pos_unf)

        # Collect results for subtrajectories (nblocks)
        v = []
        for sl in partition(self.trajectory, self.nblocks):
            grid, value = gcf_offset(f, self._discrete_tgrid, self.skip,
                                     self.trajectory.steps[sl], self._pos_unf[sl])
            v.append(value)

        # Compute variance to provide diffusion coefficient fit with weights
        self.var = numpy.std(v, axis=0) #[x if x>0 else 1e-100 for x in numpy.std(v, axis=0)]

        # Update real time grid
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]

    def analyze(self):
        # Get the time when MSD equals sigma**2
        try:
            from .helpers import feqc
            self.results['diffusive time tau_D'] = feqc(self.grid, self.value, self.sigma**2)[0]
        except:
            self.results['diffusive time tau_D'] = None

        where = (numpy.array(self.value) > self.sigma**2) * \
            (numpy.array(self.value) < self.sigma_max**2)

        if list(where).count(True) < 2:
            log.warn('could not fit MSD: not enough data above sigma')
            return

        try:
            from .helpers import linear_fit
        except ImportError:
            log.warn('could not fit MSD: missing fitting modules')
            return

        diffusion = linear_fit(numpy.array(self.grid)[where],
                               numpy.array(self.value)[where])
        ndim = self.trajectory.read(0).number_of_dimensions
        self.results['diffusion coefficient D'] = diffusion[0] / (2*ndim)
