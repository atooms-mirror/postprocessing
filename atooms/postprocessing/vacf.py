# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

""" """

import numpy

from .correlation import Correlation, gcf_offset
from .helpers import setup_t_grid

__all__ = ['VelocityAutocorrelation']


class VelocityAutocorrelation(Correlation):

    def __init__(self, trajectory, tgrid):
        Correlation.__init__(self, trajectory, tgrid, 't', 'vacf',
                             'velocity autocorrelation Z(t)', ['vel'])
        self._discrete_tgrid = setup_t_grid(trajectory, tgrid)

    def _compute(self):
        def f(x, y):
            return numpy.sum(x*y) / float(x.shape[0])
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid,
                                           self.trajectory.block_size,
                                           self.trajectory.steps,
                                           self._vel)
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]