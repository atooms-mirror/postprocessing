# This file is part of atooms
# Copyright 2010-2018, Daniele Coslovich

"""Base correlation function."""

import os
import logging
from collections import defaultdict

import numpy
from atooms.trajectory import Trajectory
from atooms.trajectory.decorators import Unfolded
from atooms.core.utils import Timer
try:
    from medepy.metadata import dump as _dump
except ImportError:
    from .helpers import _dump

from . import core
from .helpers import adjust_skip
from .progress import progress

__all__ = ['acf', 'gcf', 'gcf_offset', 'Correlation']

_log = logging.getLogger(__name__)
#pp_output_path = '{trajectory.filename}.pp.{symbol}.{tag}'
#pp_trajectory_format = None


def acf(grid, skip, t, x):
    """
    Auto correlation function.

    Calculate the correlation between time t(i0) and t(i0+i)
    for all possible pairs (i0,i) provided by grid.
    """
    cf = defaultdict(float)
    cnt = defaultdict(int)
    xave = numpy.average(x)
    for i in grid:
        for i0 in range(0, len(x)-i, skip):
            # Get the actual time difference
            dt = t[i0+i] - t[i0]
            cf[dt] += (x[i0+i]-xave) * (x[i0]-xave)
            cnt[dt] += 1

    # Return the ACF with the time differences sorted
    dt = sorted(cf.keys())
    return dt, [cf[ti] / cnt[ti] for ti in dt], cnt


def gcf(f, grid, skip, t, x):
    """
    Generalized correlation function.

    Pass a function to apply to the data.
    Exemple: mean square displacement.
    """
    # Calculate the correlation between time t(i0) and t(i0+i)
    # for all possible pairs (i0,i) provided by grid
    cf = defaultdict(float)
    cnt = defaultdict(int)
    for i in grid:
        # Note: len(x) gives x.shape[0]
        for i0 in progress(range(0, len(x)-i-1, skip)):
            # Get the actual time difference
            dt = t[i0+i] - t[i0]
            cf[dt] += f(x[i0+i], x[i0])
            cnt[dt] += 1

    # Return the ACF with the time differences sorted
    dt = sorted(cf.keys())
    return dt, [cf[ti] / cnt[ti] for ti in dt], [cnt[ti] for ti in dt]


def gcf_offset(f, grid, skip, t, x, mask=None):
    """
    Generalized correlation function

    Pass a function `f` to apply to the data `x`.
    Optionally, filter the entries at time `t[i0]` according to `mask[i0]`.

    Exemple: mean square displacement.
    """
    # Calculate the correlation between time t(i0) and t(i0+i)
    # for all possible pairs (i0,i) provided by grid
    if mask is None:
        cf = defaultdict(float)
        cnt = defaultdict(int)
        # Standard calculation
        for off, i in progress(grid, total=len(grid)):
            for i0 in range(off, len(x)-i-skip, skip):
                # Get the actual time difference
                dt = t[i0+i] - t[i0]
                cf[dt] += f(x[i0+i], x[i0])
                cnt[dt] += 1

        # Return the ACF with the time differences sorted
        dt = sorted(cf.keys())
        return dt, [cf[ti] / cnt[ti] for ti in dt]

    else:
        cf = defaultdict(float)
        cnt = defaultdict(list)
        # Filter data at time t_0 according to boolean mask
        for off, i in grid:
            for i0 in range(off, len(x)-i-skip, skip):
                # Get the actual time difference
                dt = t[i0+i] - t[i0]
                cf[dt] += f(x[i0+i][mask[i0]], x[i0][mask[i0]])
                cnt[dt].append(1)  #len(mask[i0]))

        # Return the ACF with the time differences sorted
        dt = sorted(cf.keys())
        return dt, [cf[ti] / sum([cnt[ti] for ti in dt])]


class Correlation(object):
    """
    Base class for correlation functions.

    A correlation function A(x) is defined over a grid of real
    entries {x_i} given by the list `grid`. To each entry of the grid,
    the correlation function has a corresponding value A_i=A(x_i). The
    latter values are stored in the `value` list.

    Correlation functions depending on multiple variables must provide
    a list of grids, one for each variable. The order is one specified
    by in the `symbol` variable, see below.

    Subclasses must provide a symbolic expression of the correlation
    function through the `symbol` string variable. The following
    convention is used: if a correlation function A depends on
    variables x and y, then `symbol='A(x,y)'`.

    The interface of subclasses should enforce the following
    convention for grid and samples variables passed to __init__()
    methods: if the correlation function variable is `x`, e.g. the
    correlation function is `A(x)`, then the associated grid variable
    is `xgrid` and the number of samples in that grid is `tsamples`.

    The `phasespace` variable allow subclasses to access a list of 2d
    numpy arrays with particles coordinates via the following private
    variables:

    - if `nbodies` is 1: `self._pos` for positions, `self._vel` for
    velocities, `self._unf_pos` for PBC-unfolded positions

    - if `nbodies` is 2, an additional suffix that take values 0 or 1
    is added to distinguish the two sets of particles,
    e.g. `self._pos_0` and `self._pos_1`
    """

    nbodies = 1
    """
    Class variable that controls the number of bodies, i.e. particles,
    associated to the correlation function. This variable controls the
    internal arrays used to compute the correlation, see `phasespace`.
    """
    static = False
    """
    Turn this to `True` is the correlation function is static,
    i.e. not time-dependent. This may enable some optimizations.
    """
    symbol = ''
    """Example: fskt"""
    short_name = ''
    """Example: F_s(k,t)"""
    long_name = ''
    """Example: Self intermediate scattering function"""
    phasespace = ['pos', 'pos-unf', 'vel']
    """
    List of strings or string among ['pos', 'pos-unf', 'vel']. It
    indicates which variables should be read from the trajectory file.
    They will be available as self._pos, self._pos_unf, self._vel.
    """

    def __init__(self, trj, grid, output_path=None, norigins=None):
        # Accept a trajectory-like instance or a path to a trajectory
        if isinstance(trj, str):
            self.trajectory = Trajectory(trj, mode='r', fmt=core.pp_trajectory_format)
        else:
            self.trajectory = trj
        self.grid = grid
        self.value = []
        self.analysis = {}
        self.comments = None  # can be modified by user at run time
        self.tag = ''
        self.tag_description = 'the whole system'
        self.output_path = output_path if output_path is not None else core.pp_output_path
        self.skip = adjust_skip(self.trajectory, norigins)

        # Callbacks
        self._cbk = []
        self._cbk_args = []
        self._cbk_kwargs = []

        # Lists for one body correlations
        self._pos = []
        self._vel = []
        self._pos_unf = []

        # Lists for two-body correlations
        self._pos_0, self._pos_1 = [], []
        self._vel_0, self._vel_1 = [], []
        self._pos_unf_0, self._pos_unf_1 = [], []

    def __str__(self):
        return '{} at <{}>'.format(self.long_name, id(self))

    def add_filter(self, cbk, *args, **kwargs):
        """Add filter callback `cbk` along with positional and keyword arguments"""
        if len(self._cbk) > self.nbodies:
            raise ValueError('number of filters cannot exceed n. of bodies')
        self._cbk.append(cbk)
        self._cbk_args.append(args)
        self._cbk_kwargs.append(kwargs)

    def need_update(self):
        """Check if the trajectory file is newer than the output file"""
        need = True
        if os.path.exists(self._output_file) and self._output_file != '/dev/stdout':
            if os.path.getmtime(self.trajectory.filename) < \
               os.path.getmtime(self._output_file):
                need = False
        return need

    def _setup_arrays(self):
        """
        Dump positions and/or velocities at different time frames as a
        list of numpy array.
        """
        # Ensure phasespace is a list.
        # It may not be a class variable anymore after this
        if not isinstance(self.phasespace, list) and \
           not isinstance(self.phasespace, tuple):
            self.phasespace = [self.phasespace]

        # Setup arrays
        if self.nbodies == 1:
            self._setup_arrays_onebody()
        elif self.nbodies == 2:
            self._setup_arrays_twobody()

    def _setup_arrays_onebody(self):
        """Setup list of numpy arrays for one-body correlations."""
        if 'pos' in self.phasespace or 'vel' in self.phasespace:
            for s in progress(self.trajectory):
                # Apply filter if there is one
                if len(self._cbk) > 0:
                    s = self._cbk[0](s, *self._cbk_args[0], **self._cbk_kwargs[0])
                if 'pos' in self.phasespace:
                    self._pos.append(s.dump('pos'))
                if 'vel' in self.phasespace:
                    self._vel.append(s.dump('vel'))

        # Dump unfolded positions if requested
        if 'pos-unf' in self.phasespace:
            for s in progress(Unfolded(self.trajectory, fixed_cm=True)):
                # Apply filter if there is one
                if len(self._cbk) > 0:
                    s = self._cbk[0](s, *self._cbk_args[0], **self._cbk_kwargs[0])
                self._pos_unf.append(s.dump('pos'))

    def _setup_arrays_twobody(self):
        """Setup list of numpy arrays for two-body correlations."""
        if len(self._cbk) <= 1:
            self._setup_arrays_onebody()
            self._pos_0 = self._pos
            self._pos_1 = self._pos
            return

        if 'pos' in self.phasespace or 'vel' in self.phasespace:
            for s in progress(self.trajectory):
                s0 = self._cbk[0](s, *self._cbk_args[0], **self._cbk_kwargs[0])
                s1 = self._cbk[1](s, *self._cbk_args[1], **self._cbk_kwargs[1])
                if 'pos' in self.phasespace:
                    self._pos_0.append(s0.dump('pos'))
                    self._pos_1.append(s1.dump('pos'))
                if 'vel' in self.phasespace:
                    self._vel_0.append(s0.dump('vel'))
                    self._vel_1.append(s1.dump('vel'))

        # Dump unfolded positions if requested
        if 'pos-unf' in self.phasespace:
            for s in progress(Unfolded(self.trajectory)):
                s0 = self._cbk[0](s, *self._cbk_args[0], **self._cbk_kwargs[0])
                s1 = self._cbk[1](s, *self._cbk_args[1], **self._cbk_kwargs[1])
                self._pos_unf_0.append(s0.dump('pos'))
                self._pos_unf_1.append(s1.dump('pos'))

    def compute(self):
        """
        Compute the correlation function.

        It wraps the _compute() method implemented by subclasses.
        This method also returns the tuple `self.grid`, `self.value`.
        """
        _log.info('setup arrays for %s', self.tag_description)
        t = [Timer(), Timer()]
        t[0].start()
        self._setup_arrays()
        t[0].stop()

        _log.info('computing %s for %s', self.long_name, self.tag_description)
        _log.info('using %s time origins out of %s',
                  len(range(0, len(self.trajectory), self.skip)),
                  len(self.trajectory))
        t[1].start()
        self._compute()
        t[1].stop()

        _log.info('output file %s', self._output_file)
        _log.info('done %s for %s in %.1f sec [setup:%.0f%%, compute: %.0f%%]',
                  self.long_name,
                  self.tag_description, t[0].wall_time + t[1].wall_time,
                  t[0].wall_time / (t[0].wall_time + t[1].wall_time) * 100,
                  t[1].wall_time / (t[0].wall_time + t[1].wall_time) * 100)
        _log.info('')

        return self.grid, self.value

    def _compute(self):
        """Subclasses must implement this"""
        pass

    def analyze(self):
        """
        Subclasses may implement this and store the results in the
        self.analysis dictonary
        """
        pass

    @property
    def _output_file(self):
        """Returns path of output file"""
        # Interpolate the output path string
        if self.output_path is None:
            filename = None
        else:
            filename = self.output_path.format(symbol=self.symbol,
                                               short_name=self.short_name,
                                               long_name=self.long_name.replace(' ', '_'),
                                               tag=self.tag,
                                               tag_description=self.tag_description.replace(' ', '_'),
                                               trajectory=self.trajectory)
            # Strip unpleasant punctuation
            for punct in ['.', '_', '-']:
                filename = filename.replace(punct * 2, punct)
                filename = filename.strip(punct)
        return filename

    def read(self):
        """Read correlation function from existing file"""
        with open(self._output_file, 'r') as inp:
            x = numpy.loadtxt(inp, unpack=True)
            if len(x) == 3:
                raise ValueError('cannot read 3-columns files yet')
            elif len(x) == 2:
                self.grid, self.value = x
            else:
                self.grid, self.value = x[0: 2]
                _log.warn("Ignoring some columns in %s", self._output_file)

    @property
    def grid_name(self):
        """
        Return the name of the grid variables

        Example:
        -------
        If `self.name` is `F_s(k,t)`, the function returns `['k', 't']`
        """
        variables = self.short_name.split('(')[1][:-1]
        return variables.split(',')

    def write(self):
        """Write the correlation function and the analysis data to file"""
        def is_iterable(maybe_iterable):
            try:
                iter(maybe_iterable)
            except TypeError:
                return False
            else:
                return True

        # Pack grid and value into arrays to dump
        if is_iterable(self.grid[0]) and len(self.grid) == 2:
            x = numpy.array(self.grid[0]).repeat(len(self.value[0]))
            y = numpy.array(self.grid[1] * len(self.grid[0]))
            z = numpy.array(self.value).flatten()
            dump = numpy.transpose(numpy.array([x, y, z]))
        else:
            dump = numpy.transpose(numpy.array([self.grid, self.value]))

        # Comment line
        # Extract variables from parenthesis in symbol
        variables = self.short_name.split('(')[1][:-1]
        variables = variables.split(',')
        columns = variables + [self.symbol]
        if len(self.tag_description) > 0:
            conj = 'of'
        else:
            conj = ''
        comments = _dump(title='%s %s %s %s' % (self.long_name, self.short_name, conj, self.tag_description),
                         columns=columns,
                         command='atooms-pp', version=core.__version__,
                         description=None, note=None,
                         parents=self.trajectory.filename,
                         inline=False)
        if self.comments is not None:
            comments += self.comments

        # Results of analysis
        analysis = ""
        for x, f in self.analysis.items():
            if f is not None:
                analysis += '# %s: %s\n' % (x, f)

        # Put it all together
        with open(self._output_file, 'w') as fh:
            fh.write(comments)
            if len(analysis) > 0:
                fh.write(analysis)
            numpy.savetxt(fh, dump, fmt="%g")
            fh.flush()

    def do(self, update=False):
        """
        Do the full template pattern: compute, analyze and write the
        correlation function.
        """
        if update and not self.need_update():
            self.read()
            return

        self.compute()

        try:
            self.analyze()
        except ImportError as e:
            _log.warn('Could not analyze due to missing modules, continuing...')
            _log.warn(e.message)

        self.write()

    def __call__(self):
        self.do()
