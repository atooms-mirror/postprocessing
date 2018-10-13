# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Free volume spectral density."""

import numpy
from atooms.trajectory import Trajectory
from atooms.trajectory.utils import is_cell_variable

from .fourierspace import FourierSpaceCorrelation, expo_sphere

__all__ = ['SpectralDensity']


class SpectralDensity(FourierSpaceCorrelation):
    """
    Free volume spectral density.

    From Zachary, Jiao, Torquato PRL 106, 178001 (2011).

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    symbol = 'ik'
    short_name = 'I(k)'
    long_name = 'spectral density'
    phasespace = 'pos'

    def __init__(self, trajectory, trajectory_radius, kgrid=None,
                 norigins=-1, nk=20, dk=0.1, kmin=-1.0, kmax=15.0,
                 ksamples=30):
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, norigins, nk,
                                         dk, kmin, kmax, ksamples)
        self._is_cell_variable = None
        # TODO: check step consistency 06.09.2017
        with Trajectory(trajectory_radius) as th:
            self._radius = [s.dump('particle.radius') for s in th]

    def _compute(self):
        nsteps = len(self._pos)
        # Setup k vectors and tabulate rho
        kgrid, selection = self.kgrid, self.selection
        kmax = max(self.kvector.keys()) + self.dk
        cnt = [0 for k in kgrid]
        # Note: actually rho_av is not calculated because it is negligible
        rho_av = [complex(0., 0.) for k in kgrid]
        rho2_av = [complex(0., 0.) for k in kgrid]
        cell_variable = is_cell_variable(self.trajectory)
        for i in range(0, nsteps, self.skip):
            # If cell changes we have to update
            if cell_variable:
                self._setup(i)
                kgrid, selection = self._decimate_k()
                kmax = max(self.kvector.keys()) + self.dk

            expo = expo_sphere(self.k0, kmax, self._pos[i])
            for kk, knorm in enumerate(kgrid):
                for k in selection[kk]:
                    ik = self.kvector[knorm][k]
                    Ri = self._radius[i]
                    mk = 4 * numpy.pi / knorm**3 * (numpy.sin(knorm*Ri) - (knorm*Ri) * numpy.cos(knorm*Ri))
                    rho = numpy.sum(mk * expo[..., 0, ik[0]] * expo[..., 1, ik[1]] * expo[..., 2, ik[2]])
                    rho2_av[kk] += (rho * rho.conjugate())
                    cnt[kk] += 1

        # Normalization.
        volume = numpy.average([s.cell.volume for s in self.trajectory])
        self.grid = kgrid
        self.value = [(rho2_av[kk] / cnt[kk] - rho_av[kk]*rho_av[kk].conjugate() / cnt[kk]**2).real / volume
                      for kk in range(len(self.grid))]
        self.value_nonorm = [rho2_av[kk].real / cnt[kk]
                             for kk in range(len(self.grid))]
