# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Real space correlation functions."""

import numpy
import math
import logging

from .helpers import linear_grid, logx_grid
from .correlation import Correlation, gcf_offset
from .helpers import adjust_skip, setup_t_grid

log = logging.getLogger(__name__)


# Kernels

def partition(inp, nbl):
    nel = len(inp) / nbl
    a = []
    for i in range(nbl):
        a.append(slice(i * nel, (i+1) * nel))
    return a

def gr_kernel(x, y, L):
    # r is an array of array distances
    r = x-y
    r = r - numpy.rint(r/L) * L
    return numpy.sqrt(numpy.sum(r**2, axis=1))

def gr_kernel_square(x, y, L):
    """Return square distances."""
    # r is an array of array distances
    r = x-y
    r = r - numpy.rint(r/L) * L
    return numpy.sum(r**2, axis=1)

def pairs(f, x, y, L):
    """Apply function f to all pairs in x[i] and y[j] and return a list of
    f values.
    """
    fxy = []
    for i in xrange(len(y)):
        fxy.append(f(x[:], y[i], L))
    return fxy

def pairs_numpy(f, x, y, L):
    """Apply function f to all pairs in x[i] and y[j] and return a numpy
    array of f values.
    """
    fxy = numpy.ndarray((len(y), len(x)))
    for i in xrange(fxy.shape[0]):
        fxy[i, :] = f(x[:], y[i], L)
    return fxy

def pairs_newton(f, x, y, L):
    """Apply function f to all pairs in x[i] and y[j] and return a list of
    f values.
    """
    fxy = []
    for i in xrange(len(y)-1):
        for value in f(x[i+1:], y[i], L):
            fxy.append(value)
    return fxy

def pairs_newton_hist(f, x, y, L, bins):
    """Apply function f to all pairs in x[i] and y[j] and update the
    |hist| histogram using the |bins| bin edges.
    """
    hist, bins = numpy.histogram([], bins)
    # Do the calculation in batches to optimize
    bl = max(1, int(100 * 1000.0 / len(y)))
    for ib in range(0, len(y)-1, bl):
        fxy = []
        # batch must never exceed len(y)-1
        for i in xrange(ib, min(ib+bl, len(y)-1)):
            for value in f(x[i+1:], y[i], L):
                fxy.append(value)
        hist_tmp, bins = numpy.histogram(fxy, bins)
        hist += hist_tmp
    return hist

def pairs_hist(f, x, y, L, bins):
    """Apply function f to all pairs in x[i] and y[j] and update the
    |hist| histogram using the |bins| bin edges.
    """
    hist, bins = numpy.histogram([], bins)
    for i in xrange(len(y)):
        fxy = f(x[:], y[i], L)
        hist_tmp, bins = numpy.histogram(fxy, bins)
        hist += hist_tmp
    return hist

def square_displacement(x, y, L=None):
    """Return array of square distances (no pbc)."""
    return numpy.sum((x-y)**2, axis=1)

def mean_square_displacement(x, y):
    return numpy.sum((x-y)**2) / float(x.shape[0])

def non_gaussian_parameter(x, y):
    if x is y:
        return 0.0
    dx2 = (x-y)**2
    dr2 = numpy.sum(dx2) / float(x.shape[0])
    dr4 = numpy.sum(numpy.sum(dx2, axis=1)**2) / float(x.shape[0])
    return 3*dr4 / (5*dr2**2) - 1

def collective_overlap(r0, r1, side, a_square):
    rij = pairs_numpy(square_displacement, r0, r1, side)
    return (rij.flatten() < a_square).sum()

def self_overlap(r0, r1, side, a_square):
    rij = square_displacement(r0, r1)
    return rij.flatten() < a_square


# Correlation classes

class CollectiveOverlap(Correlation):

    # TODO: why dont we use PBCed distance here?!

    def __init__(self, trajectory, grid=None, nsamples=60, a=0.3,
                 norigins=-1):
        Correlation.__init__(self, trajectory, grid, 't', 'qt',
                             'collective overlap Q(t)', 'pos')
        self.a_square = a**2
        self.skip = adjust_skip(self.trajectory, norigins)
        if grid is None:
            self.grid = logx_grid(0.0, trajectory.total_time * 0.75, nsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid)

    def _compute(self):
        side = self.trajectory.read(0).cell.side
        def f(x, y):
            return collective_overlap(x, y, side, self.a_square).sum() / float(x.shape[0])        
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid,
                                           self.skip, self.trajectory.steps, self._pos)
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]


class SelfOverlap(Correlation):

    def __init__(self, trajectory, grid=None, norigins=-1, a=0.3,
                 nsamples=60):
        Correlation.__init__(self, trajectory, grid, 't', 'qst',
                             'self overlap Q_s(t)', 'pos-unf')
        if not self._need_update:
            return
        if grid is None:
            self.grid = logx_grid(0.0, trajectory.total_time * 0.75, nsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid)
        self.skip = adjust_skip(self.trajectory, norigins)
        self.a_square = a**2

    def _compute(self):
        side = self.trajectory.read(0).cell.side
        def f(x, y):
            return self_overlap(x, y, side, self.a_square).sum() / float(x.shape[0])
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid,
                                           self.skip, self.trajectory.steps, self._pos_unf)
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]

    def analyze(self):
        try:
            from .helpers import feqc
            self.results['tau'] = feqc(self.grid, self.value, 1/numpy.exp(1.0))[0]
        except:
            self.results['tau'] = None


class Chi4SelfOverlap(Correlation):

    # TODO: refactor correlation init via class variables??
    def __init__(self, trajectory, grid=None, norigins=-1, a=0.3,
                 nsamples=60):
        Correlation.__init__(self, trajectory, grid, 't', 'chi4qs',
                             'dynamic scusceptibility of self overlap chi_4(t)', 'pos-unf')
        if not self._need_update:
            return
        if grid is None:
            self.grid = logx_grid(0.0, trajectory.total_time * 0.75, nsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid)
        self.skip = adjust_skip(self.trajectory, norigins)
        self.a_square = a**2
        self.average = Correlation(trajectory, self.grid, 't', 'qsu',
                                   'Average of self overlap not normalized')
        self.variance = Correlation(trajectory, self.grid, 't', 'qs2u',
                                    'Variance self overlap not normalized')

    def _compute(self):
        # TODO: write general susceptibility
        # At this stage, we must copy over the tags
        self.average.tag, self.variance.tag = self.tag, self.tag
        side = self.trajectory.read(0).cell.side
        def f(x, y):
            return self_overlap(x, y, side, self.a_square).sum()

        self.grid = []
        for off, i in self._discrete_tgrid:
            A, A2, cnt = 0.0, 0.0, 0
            for i0 in xrange(off, len(self._pos_unf)-i-self.skip, self.skip):
                w = f(self._pos_unf[i0], self._pos_unf[i0+i])
                A2 += w**2
                A += w
                cnt += 1
            dt = self.trajectory.steps[off+i] - self.trajectory.steps[off]
            A_av = A/cnt
            A2_av = A2/cnt
            self.grid.append(dt * self.trajectory.timestep)
            self.value.append((A2_av - A_av**2) / self._pos_unf[0].shape[0])
            self.average.value.append(A_av)
            self.variance.value.append(A2_av)
        self.average.grid, self.variance.grid = self.grid, self.grid

    def write(self):
        # We subclass this to also write down qsu and qsu2
        super(Chi4SelfOverlap, self).write()
        self.average.write()
        self.variance.write()

    def analyze(self):
        from .helpers import ifabsmm
        try:
            self.results['tau_star'], self.results['chi4_star'] = ifabsmm(self.grid, self.value)[1]
        except ZeroDivisionError:
            print '# warning : could not find maximum'
            pass


class OverlapDistribution(Correlation):

    # TODO: OverlapDistribution is not really a Correlation but more an Analysis kind of object
    # TODO: filter matrix particles from overlap !

    def __init__(self, trajectory, grid, skip=1, a=0.3):
        Correlation.__init__(self, trajectory, grid, 'q', 'P_q',
                             'overlap distribution P(q)', 'pos')
        self.skip = skip
        self.a = a
        self.a_square = a**2

    def _compute(self):
        side = self.trajectory.read(0).cell.side
        N = len(self.trajectory.read(0).particle)
        #N_m = len(self.trajectory.read(0).matrix)
        data = []
        self_threshold = 10.0 #/ N
        for i0 in range(0, len(self._pos), self.skip):
            for i in range(i0, len(self._pos), self.skip):
                qs = self_overlap(self._pos[i0], self._pos[i],
                                  side, self.a_square).sum() #/ float(N)

                # We apply at every time step because a replica might have come back
                # to the initial state without having travelled around long enough
                if qs < self_threshold:
                    q = collective_overlap(self._pos[i0], self._pos[i],
                                           side, self.a_square) #/ float(N)
                    data.append(q)

        if len(data) == 0:
            print 'No data matched overlap criteria'

        # TODO: normalization and binning could be done later, in a script
        #norm_data = numpy.array(data) / (float(N) * float(N_m)/(N+N_m))
        norm_data = numpy.array(data) / float(N) #- float(N_m)/N
        self.value, self.grid = numpy.histogram(norm_data, bins=self.grid, normed=True)
        self.grid = [(x1 + x2)/2.0 for x1, x2 in zip(self.grid[0:], self.grid[1:])]


def block_matrix(n, weight, max_size):
    """ Determine optimal partitioning of a nxn matrix into blocks of maximal area max_size """
    nbl = max(1, int(n / float(max_size)**0.5)) * weight
    nbl = min(nbl, n)
    grid = linear_grid(0, n, nbl)
    # make sure the last entry is n-1
    if grid[-1] != n-1:
        grid += [n-1]
    return grid


class RadialDistributionFunction(Correlation):

    nbodies = 2

    def __init__(self, trajectory, grid=None, norigins=-1, dr=0.04):
        Correlation.__init__(self, trajectory, grid, 'r', 'gr',
                             'radial distribution function g(r)', 'pos')
        self.skip = adjust_skip(trajectory, norigins)
        self.side = self.trajectory.read(0).cell.side
        if grid is not None:
            # Reconstruct bounds of grid for numpy histogram
            self.grid = []
            for i in range(len(grid)):
                self.grid.append(grid[i] - (grid[1]-grid[0])/2)
            self.grid.append(grid[-1] + (grid[1]-grid[0])/2)
        else:
            self.grid = linear_grid(0.0, self.side[0]/2.0, dr)

    def _compute(self):
        ncfg = len(self.trajectory)
        if self.trajectory.grandcanonical:
            N_0 = numpy.average([len(x) for x in self._pos_0])
            N_1 = numpy.average([len(x) for x in self._pos_1])
        else:
            N_0, N_1 = len(self._pos_0[0]), len(self._pos_1[0])

        gr_all = []
        _, r = numpy.histogram([], bins=self.grid)
        for i in range(0, ncfg, self.skip):
            if self._pos_0 is self._pos_1:
                gr = pairs_newton_hist(gr_kernel, self._pos_0[i], self._pos_1[i],
                                       self.side, r)
            else:
                gr = pairs_hist(gr_kernel, self._pos_0[i], self._pos_1[i],
                                self.side, r)
            gr_all.append(gr)

        # Normalization
        vol = 4 * math.pi / 3.0 * (r[1:]**3-r[:-1]**3)
        rho = N_0 / self.side.prod()
        if self._pos_0 is self._pos_1:
            norm = rho * vol * N_0 * 0.5  # use Newton III
        else:
            norm = rho * vol * N_1
        gr = numpy.average(gr_all, axis=0)
        self.grid = (r[:-1]+r[1:]) / 2.0
        self.value = gr / norm


class MeanSquareDisplacement(Correlation):

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
            self.results['tau_D'] = feqc(self.grid, self.value, self.sigma**2)[0]
        except:
            self.results['tau_D'] = None

        where = (numpy.array(self.value) > self.sigma**2) * \
            (numpy.array(self.value) < self.sigma_max**2)

        if list(where).count(True) < 2:
            log.warn('could not fit MSD: not enough data above sigma')
            return

        try:
            from .helpers import linear_fit
        except ImportError as e:
            log.warn('could not fit MSD: missing fitting modules')
            return

        diffusion = linear_fit(numpy.array(self.grid)[where],
                               numpy.array(self.value)[where])
        ndim = self.trajectory.read(0).number_of_dimensions
        self.results['diffusion coefficient D'] = diffusion[0] / (2*ndim)


class NonGaussianParameter(Correlation):

    def __init__(self, trajectory, tgrid=None, norigins=50, nsamples=30):
        Correlation.__init__(self, trajectory, tgrid, 't', 'alpha2',
                             "non-Gaussian parameter alpha_2(t)", ['pos-unf'])
        if not self._need_update:
            return
        if self.grid is None:
            self.grid = linear_grid(0.0, trajectory.total_time * 0.75, nsamples)
        self._discrete_tgrid = setup_t_grid(trajectory, self.grid)
        self.skip = adjust_skip(trajectory, norigins)

    def _compute(self):
        f = non_gaussian_parameter
        self.grid, self.value = gcf_offset(f, self._discrete_tgrid, self.skip,
                                           self.trajectory.steps, self._pos_unf)
        self.grid = [ti * self.trajectory.timestep for ti in self.grid]

    def analyze(self):
        try:
            from .helpers import ifabsmm
            self.results['t_star'], self.results['a2_star'] = ifabsmm(self.grid, self.value)[1]
        except:
            pass


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
