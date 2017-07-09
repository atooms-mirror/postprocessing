#!/usr/bin/env python

"""
Bond breaking probability.

Compute the probability of loosing all but `n` bonds (neighbors) after
a time t.

Input file must contain neighbors information.
"""

import sys
import datetime
import atooms.trajectory as trj
import atooms.voronoi as vor
from postprocessing.helpers import setup_t_grid, adjust_skip, logx_grid
from postprocessing.correlation import gcf_offset

def main(f, bonds=1, norigins=40):

    with vor.TrajectoryVoronoi(f) as th:

        # Put data into a big list of lists
        npart = len(th[0].particle)
        neighbors_list = []
        for frame, s in enumerate(th):
            data = []
            for i in range(npart):
                nn = set(s.particle[i].neighbors)
                data.append(nn)
            neighbors_list.append(data)

        # Kernel
        def func(x, y):
            broken = [int(len(x[i] & y[i]) > bonds) for i in range(npart)]
            return sum(broken) / float(npart)

        # Compute probability
        skip = adjust_skip(th, -1)
        grid = logx_grid(0.0, th.total_time * 0.75, norigins)
        discrete_tgrid = setup_t_grid(th, grid)
        grid, value = gcf_offset(func, discrete_tgrid, skip, th.steps,
                                 neighbors_list)

        # Dump
        with open('%s.bondbreak-%d' % (f, bonds), 'w') as fh:
            fh.write('# title: probability P_b(t) of loosing all but n=%s bonds after a time t\n' % bonds)
            fh.write('# columns: t, P_b(t)\n')
            fh.write('# created: %s\n' % datetime.datetime.now())
            fh.write('# notes: postprocessing v%s\n' % 0.1)
            for x, y in zip(grid, value):
                fh.write('%g %g\n' % (x, y))

import argh
argh.dispatch_command(main)

