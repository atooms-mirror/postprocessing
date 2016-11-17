# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

#!/usr/bin/env python

import sys
import numpy
from .helpers import linear_grid, logx_grid
from atooms import trajectory
from .correlation import Correlation
from . import correlation

correlation.LOG_LEVEL = 'INFO'

class StressAutocorrelation(Correlation):
    
    def __init__(self, trajectory, tgrid):
        Correlation.__init__(self, trajectory, tgrid, 't', 'sacf', "Stress autocorrelation", ['vel'])
        self._discrete_tgrid = correlation._setup_t_grid(trajectory, tgrid)

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
        
        # for t, s in zip(self.trajectory.steps, self._stress):
        #     print t, s

        self.grid, self.value = correlation.gcf_offset(f, self._discrete_tgrid, self.trajectory.block_period, 
                                                       self.trajectory.steps, self._stress)

        self.value = [x / V for x in self.value]
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]

    def analyze(self):
        pass

def main(f):
#    t = trajectory.Trajectory('/home/coslo/reference/kalj/config.dat')
#    t = trajectory.Trajectory('/home/coslo/projects/sandwich/data/lj-fluid-nvt-ncfg1000_boost4-exp/elong1.00/nl5/rho1.0000_T2.100/n00/config.dat')
    t = trajectory.Trajectory(f)
    s = StressAutocorrelation(t, linear_grid(0.0, 1.5, 60))
#    s = StressAutocorrelation(t, logx_grid(0.0, t.time_total/10.0, 60))
    s.do()

if __name__ == '__main__':
    for f in sys.argv[1:]:
        main(f)
    
