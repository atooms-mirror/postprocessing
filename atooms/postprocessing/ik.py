# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Free volume spectral density."""

import numpy

from .helpers import adjust_skip
from .fourierspace import FourierSpaceCorrelation, expo_sphere
from atooms.trajectory.utils import is_cell_variable

__all__ = ['SpectralDensity']


class SpectralDensity(FourierSpaceCorrelation):

    """
    Free volume spectral density.

    From Zachary, Jiao, Torquato PRL 106, 178001 (2011).

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    def __init__(self, trajectory, trajectory_radius, kgrid=None,
                 norigins=-1, nk=20, dk=0.1, kmin=-1.0, kmax=15.0,
                 ksamples=30):
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, 'I(k)',
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
