#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Stress autocorrelation function."""

import numpy
from .correlation import Correlation, gcf_offset
from .helpers import setup_t_grid


class StressAutocorrelation(Correlation):

    """Stress autocorrelation function."""

    def __init__(self, trajectory, tgrid):
        Correlation.__init__(self, trajectory, tgrid, 'S(t)', 'sacf',
                             'stress autocorrelation', ['vel'])
        self._discrete_tgrid = setup_t_grid(trajectory, tgrid)

    def _get_stress(self):
        ndims = 3
        p = self.trajectory.read(0).particle
        V = self.trajectory.read(0).cell.volume
        mass = numpy.array([pi.mass for pi in p])
        self._stress = []
        for i in self.trajectory.samples:
            s = self.trajectory.read(i).interaction.total_stress
            l = 0
            # Must be recreated otherwise we just append references to it
            slk = numpy.zeros(ndims)
            for j in range(ndims):
                for k in range(j+1,ndims):
                    slk[l] = s[j,k] + numpy.sum(mass[:] * self._vel[i][:,j] * self._vel[i][:,k])
#                    slk[l] = numpy.sum(mass[:] * self._vel[i][:,j] * self._vel[i][:,k])
                    l += 1
            self._stress.append(slk)

    def _compute(self):
        def f(x, y):
            return numpy.sum(x*y) / float(x.shape[0])

        self._get_stress()
        V = self.trajectory.read(0).cell.volume
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid, self.trajectory.block_size,
                                           self.trajectory.steps, self._stress)
        self.value = [x / V for x in self.value]
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]

    def analyze(self):
        pass
