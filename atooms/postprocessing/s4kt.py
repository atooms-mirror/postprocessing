# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Four-point dynamic structure factor."""

import numpy

from .fourierspace import FourierSpaceCorrelation, expo_sphere
from .helpers import adjust_skip, setup_t_grid
from .qt import self_overlap

__all__ = ['S4ktOverlap']


class S4ktOverlap(FourierSpaceCorrelation):

    """
    Four-point dynamic structure factor from time-dependent self overlap.

    See the documentation of the `FourierSpaceCorrelation` base class
    for information on the instance variables.
    """

    # TODO: refactor a S4k base correlation that forces to implement tabulat method (e.g. overlap, Q_6, voronoi ...)
    # TODO: should we drop this instead and rely on F(k,t) with grandcanonical

    def __init__(self, trajectory, tgrid, kgrid=None, norigins=-1,
                 nk=20, dk=0.1, a=0.3, kmin=1.0, kmax=10.0, ksamples=10):
        FourierSpaceCorrelation.__init__(self, trajectory, [tgrid, kgrid], 'S_4(k,t)', 's4kt',
                                         '4-point dynamic structure factor from self overlap',
                                         ['pos', 'pos-unf'],
                                         nk, dk, kmin, kmax, ksamples)
        # Setup time grid
        self._discrete_tgrid = setup_t_grid(trajectory, tgrid)
        self.skip = adjust_skip(self.trajectory, norigins)
        self.a_square = a**2

        # Setup k vectors and tabulate rho
        # TODO: move decimate up the chain?
        self.k_sorted, self.k_selected = self._decimate_k()
        # Redefine kgrid to give exactly the average wave vectors used.
        # TODO; should we do it for in base?
        self.grid[1] = self._actual_k_grid(self.k_sorted, self.k_selected)

    def _tabulate_W(self, k_sorted, k_selected, t_off, t, skip):
        """ Tabulate W """
        nsteps = len(self._pos)
        side = self.trajectory[0].cell.side
        kmax = max(self.kvec.keys()) + self.dk
        nt = range(t_off, len(self._pos)-t, skip)
        W = {}
        for i_0, t_0 in enumerate(nt):
            expo = expo_sphere(self.k0, kmax, self._pos[t_0])
            for kk, knorm in enumerate(k_sorted):
                for i in k_selected[kk]:
                    ik = self.kvec[knorm][i]
                    if not ik in W:
                        W[ik] = numpy.ndarray(len(nt), dtype=complex)
                    W[ik][i_0] = numpy.sum(self_overlap(self._pos_unf[t_0], self._pos_unf[t_0+t], side, self.a_square) *
                                          expo[...,0,ik[0]] * expo[...,1,ik[1]] * expo[...,2,ik[2]])
        return W

    def _compute(self):
        # Make sure there is only one time in tgrid.
        # We could easily workaround it by outer looping over i
        # We do not expect to do it for many times (typically we show S_4(k,tau_alpha) vs k)
        # if len(self._discrete_tgrid) > 1:
        #     raise ValueError('There should be only one time for S4kt')
        dt = []
        self.value = []
        for off, i  in self._discrete_tgrid:

            # as for fkt
            W = self._tabulate_W(self.k_sorted, self.k_selected, off, i, self.skip)

            # Compute vriance of W
            cnt = [0 for k in self.k_sorted]
            w_av = [complex(0., 0.) for k in self.k_sorted]
            w2_av = [complex(0., 0.) for k in self.k_sorted]
            for kk, knorm in enumerate(self.k_sorted):
                for j in self.k_selected[kk]:
                    ik = self.kvec[knorm][j]
                    # Comupte |<W>|^2  and <W W^*>
                    w_av[kk] = numpy.average(W[ik])
                    w2_av[kk] = numpy.average(W[ik] * W[ik].conjugate())

            # Normalization
            npart = self._pos[0].shape[0]
            dt.append(self.trajectory.timestep * (self.trajectory.steps[off+i] - self.trajectory.steps[off]))
            #self.grid[1] = k_sorted
            self.value.append([float(w2_av[kk] - (w_av[kk]*w_av[kk].conjugate())) / npart for kk in range(len(self.grid[1]))])
        self.grid[0] = dt
