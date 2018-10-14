# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Structure factor."""

import logging
import numpy

from .progress import progress
from .fourierspace import FourierSpaceCorrelation, expo_sphere

__all__ = ['StructureFactor', 'StructureFactorOptimized']

_log = logging.getLogger(__name__)


class StructureFactor(FourierSpaceCorrelation):
    """
    Structure factor.

    If `trajectory_field` is not `None`, the field is read from the
    last column of this trajectory file, unless the `field` string is
    provided.

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    nbodies = 2
    symbol = 'sk'
    short_name = 'S(k)'
    long_name = 'structure factor'
    phasespace = ['pos']

    def __init__(self, trajectory, kgrid=None, norigins=-1, nk=20,
                 dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30,
                 trajectory_field=None, field=None):
        FourierSpaceCorrelation.__init__(self, trajectory, kgrid, norigins,
                                         nk, dk, kmin, kmax, ksamples)
        self._is_cell_variable = None
        self._field, tag = self._add_field(trajectory_field, field)
        if tag is not None:
            self.tag = tag
            self.tag_description += '%s field' % tag.replace('_', ' ')

    def _add_field(self, trajectory_field, field):
        if trajectory_field is None:
            return None, None
        else:
            # TODO: check step consistency 06.09.2017
            from atooms.trajectory import TrajectoryXYZ
            with TrajectoryXYZ(trajectory_field) as th:
                if th.steps != self.trajectory.steps:
                    raise ValueError('field and traectory are not synced (%s, %s)' % (len(th), len(self.trajectory)))
                fields = []
                # This must be a string, not a list
                unique_field = th._read_metadata(0)['columns']
                if isinstance(unique_field, list):
                    # If field is not given, get the last column
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
        kgrid, selection = self.kgrid, self.selection
        kmax = max(self.kvector.keys()) + self.dk
        cnt = [0 for k in kgrid]
        rho_av = [complex(0., 0.) for k in kgrid]
        rho2_av = [complex(0., 0.) for k in kgrid]
        variable_cell = is_cell_variable(self.trajectory)
        for i in progress(range(0, nsteps, self.skip), total=nsteps // self.skip):
            # If cell changes we have to update the wave vectors
            if variable_cell:
                self._setup(i)
                kgrid, selection = self._decimate_k()
                kmax = max(self.kvector.keys()) + self.dk

            # Tabulate exponentials
            # Note: tabulating and computing takes about the same time
            if self._pos_0[i] is self._pos_1[i]:
                # Identical species
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_0
            else:
                # Cross correlation
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_sphere(self.k0, kmax, self._pos_1[i])

            for kk, knorm in enumerate(kgrid):
                for k in selection[kk]:
                    ik = self.kvector[knorm][k]
                    # In the absence of a microscopic field, rho_av = (0, 0)
                    if not self._field:
                        if expo_0 is expo_1:
                            # Identical species
                            rho_0 = numpy.sum(expo_0[..., 0, ik[0]] *
                                              expo_0[..., 1, ik[1]] *
                                              expo_0[..., 2, ik[2]])
                            rho_1 = rho_0
                        else:
                            # Cross correlation
                            rho_0 = numpy.sum(expo_0[..., 0, ik[0]] *
                                              expo_0[..., 1, ik[1]] *
                                              expo_0[..., 2, ik[2]])
                            rho_1 = numpy.sum(expo_1[..., 0, ik[0]] *
                                              expo_1[..., 1, ik[1]] *
                                              expo_1[..., 2, ik[2]])
                    else:
                        # We have a field as a weight
                        rho_0 = numpy.sum(self._field[i] *
                                          expo_0[..., 0, ik[0]] *
                                          expo_0[..., 1, ik[1]] *
                                          expo_0[..., 2, ik[2]])
                        rho_1 = rho_0
                        rho_av[kk] += rho_0

                    rho2_av[kk] += (rho_0 * rho_1.conjugate())
                    cnt[kk] += 1

        # Normalization.
        npart_0 = sum([p.shape[0] for p in self._pos_0]) / float(len(self._pos_0))
        npart_1 = sum([p.shape[0] for p in self._pos_1]) / float(len(self._pos_1))
        self.grid = kgrid
        self.value, self.value_nonorm = [], []
        for kk in range(len(self.grid)):
            norm = float(npart_0 * npart_1)**0.5
            value = (rho2_av[kk] / cnt[kk] - rho_av[kk]*rho_av[kk].conjugate() / cnt[kk]**2).real
            self.value.append(value / norm)
            self.value_nonorm.append(value)


class StructureFactorOptimized(StructureFactor):
    """
    Optimized structure factor.

    It uses a fortran 90 extension.
    """

    nbodies = 2
    symbol = 'sk'
    short_name = 'S(k)'
    long_name = 'structure factor'
    phasespace = ['pos']

    def _compute(self):
        from atooms.trajectory.utils import is_cell_variable
        try:
            from atooms.postprocessing.fourierspace_wrap import fourierspace_module
        except ImportError:
            _log.error('f90 wrapper missing or not functioning')
            raise

        nsteps = len(self._pos_0)
        # Setup k vectors and tabulate rho
        kgrid, selection = self.kgrid, self.selection
        kmax = max(self.kvector.keys()) + self.dk
        cnt = [0 for k in kgrid]
        rho_av = [complex(0., 0.) for k in kgrid]
        rho2_av = [complex(0., 0.) for k in kgrid]
        variable_cell = is_cell_variable(self.trajectory)
        for i in range(0, nsteps, self.skip):
            # If cell changes we have to update the wave vectors
            if variable_cell:
                self._setup(i)
                kgrid, selection = self._decimate_k()
                kmax = max(self.kvector.keys()) + self.dk

            # Tabulate exponentials
            # Note: tabulating and computing takes about the same time
            if self._pos_0[i] is self._pos_1[i]:
                # Identical species
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_0
            else:
                # Cross correlation
                # TODO: cross correlation wont work
                expo_0 = expo_sphere(self.k0, kmax, self._pos_0[i])
                expo_1 = expo_sphere(self.k0, kmax, self._pos_1[i])

            for kk, knorm in enumerate(kgrid):
                ikvec = numpy.ndarray((3, len(selection[kk])), order='F', dtype=numpy.int32)
                i = 0
                for k in selection[kk]:
                    ikvec[:, i] = self.kvector[knorm][k]
                    i += 1
                rho = numpy.zeros(ikvec.shape[1], dtype=numpy.complex128)
                fourierspace_module.sk_bare(expo_0, ikvec, rho)
                rho_0 = rho
                rho_1 = rho
                rho2_av[kk] += numpy.sum(rho_0 * rho_1.conjugate())
                cnt[kk] += rho.shape[0]

        # Normalization.
        npart_0 = sum([p.shape[0] for p in self._pos_0]) / float(len(self._pos_0))
        npart_1 = sum([p.shape[0] for p in self._pos_1]) / float(len(self._pos_1))
        self.grid = kgrid
        self.value, self.value_nonorm = [], []
        for kk in range(len(self.grid)):
            norm = float(npart_0 * npart_1)**0.5
            value = (rho2_av[kk] / cnt[kk] - rho_av[kk]*rho_av[kk].conjugate() / cnt[kk]**2).real
            self.value.append(value / norm)
            self.value_nonorm.append(value)


