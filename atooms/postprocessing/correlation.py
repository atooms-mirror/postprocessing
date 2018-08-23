# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import sys
import os
import copy
import numpy
import math
import random
import warnings
import logging
from collections import defaultdict
try:
    from medepy.metadata import dump as _dump
except ImportError:
    from .helpers import _dump
from .core import __version__
from atooms.trajectory.decorators import Unfolded


log = logging.getLogger(__name__)


def acf(grid, skip, t, x):
    """Auto correlation function.
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
    return dt, [cf[t] / cnt[t] for t in dt], cnt

def gcf(f, grid, skip, t, x):
    """Generalized correlation function.
    Pass a function to apply to the data.
    Exemple: mean square displacement.
    """
    # Calculate the correlation between time t(i0) and t(i0+i)
    # for all possible pairs (i0,i) provided by grid
    acf = defaultdict(float)
    cnt = defaultdict(int)
    for i in grid:
        # Note: len(x) gives x.shape[0]
        for i0 in xrange(0, len(x)-i-1, skip):
            # Get the actual time difference
            dt = t[i0+i] - t[i0]
            acf[dt] += f(x[i0+i], x[i0])
            cnt[dt] += 1

    # Return the ACF with the time differences sorted
    dt = sorted(acf.keys())
    return dt, [acf[t] / cnt[t] for t in dt], [cnt[t] for t in dt]

def gcf_offset(f, grid, skip, t, x, mask=None):
    """
    Generalized correlation function

    Pass a function `f` to apply to the data `x`.

    Optionnally, filter the entries at time `t[i0]` according to `mask[i0]`.

    Exemple: mean square displacement.
    """
    # Calculate the correlation between time t(i0) and t(i0+i)
    # for all possible pairs (i0,i) provided by grid
    if mask is None:
        acf = defaultdict(float)
        cnt = defaultdict(int)
        # Standard calculation
        for off, i in grid:
            for i0 in xrange(off, len(x)-i-skip, skip):
                # Get the actual time difference
                dt = t[i0+i] - t[i0]
                acf[dt] += f(x[i0+i], x[i0])
                cnt[dt] += 1

        # Return the ACF with the time differences sorted
        dt = sorted(acf.keys())
        return dt, [acf[t] / cnt[t] for t in dt] #, [cnt[t] for t in dt]

    else:
        acf = defaultdict(float)
        cnt = defaultdict(list)
        # Filter data at time t_0 according to boolean mask
        for off, i in grid:
            for i0 in xrange(off, len(x)-i-skip, skip):
                # Get the actual time difference
                dt = t[i0+i] - t[i0]
                acf[dt] += f(x[i0+i][mask[i0]], x[i0][mask[i0]])
                cnt[dt].append(1)  #len(mask[i0]))

        # Return the ACF with the time differences sorted
        dt = sorted(acf.keys())
        return dt, [acf[t] / sum([cnt[t] for t in dt])] #, [cnt[t] for t in dt]


UPDATE = False
OUTPUT_PATH = '{trajectory.filename}.pp.{short_name}.{tag}'

class Correlation(object):

    nbodies = 1

    def __init__(self, trj, grid, variables='', short_name='',
                 description='', phasespace=None,
                 output_path=None):
        # TODO: we could force trajectory cast if a string is passed
        # self.variables = ('k', 't')
        # self.short_name = 'fskt'
        # self.description = 'Self intermediate scattering function F_s(k,t)'
        # self.tag_description = 'of particles A'  # 'of radius field'
        # self.tag = 'A'
        self.trajectory = trj
        self.grid = grid
        self.variables = variables
        self.short_name = short_name
        self.description = description
        self.results = {}
        self.output = None
        if phasespace is None:
            self._phasespace = []
        else:
            self._phasespace = phasespace
        self.output_path = output_path if output_path is not None else OUTPUT_PATH
        self.tag = ''
        self.tag_description = ''
        self.comments = None # can be modified by user at run time
        if isinstance(phasespace, str):
            self._phasespace = [phasespace]
        self.value = []
        self.cbk = []
        self.cbk_args = []
        self.cbk_kwargs = []
        log.debug('%s for %s' % (self.description, self.trajectory.filename))

        # If update mode is on, we will only do the calculation if the trajectory
        # file is newer than any of the provided files
        self._need_update = True
        if UPDATE:
            if os.path.exists(self._output_file) and self._output_file == '/dev/stdout':
                if os.path.getmtime(self.trajectory.filename) < \
                   os.path.getmtime(self._output_file):
                    self._need_update = False
                    # # TODO: to optimize avoid reading correlation objects unless we explicitly pass something to __init__
                    self.read()

    def __str__(self):
        return '%s' % (self.description, )

    def add_filter(self, cbk, *args, **kwargs):
        if len(self.cbk) > self.nbodies:
            raise ValueError('number of filters cannot exceed n. of bodies')
        self.cbk.append(cbk)
        self.cbk_args.append(args)
        self.cbk_kwargs.append(kwargs)

    def _setup_arrays(self):
        """Dump positions and/or velocities in numpy array"""
        if self.nbodies == 1:
            self._setup_arrays_onebody()
        elif self.nbodies == 2:
            self._setup_arrays_twobody()

    def _setup_arrays_onebody(self):
        self._pos = []
        self._vel = []
        if 'pos' in self._phasespace or 'vel' in self._phasespace:
            for s in self.trajectory:
                # Apply filter if there is one
                if len(self.cbk) > 0:
                    s = self.cbk[0](s, *self.cbk_args[0], **self.cbk_kwargs[0])
                if 'pos' in self._phasespace:
                    self._pos.append(s.dump('pos'))
                if 'vel' in self._phasespace:
                    self._vel.append(s.dump('vel'))

        # Dump unfolded positions if requested
        self._pos_unf = []
        if 'pos-unf' in self._phasespace:
            for s in Unfolded(self.trajectory, fixed_cm=True):
                # Apply filter if there is one
                if len(self.cbk) > 0:
                    s = self.cbk[0](s, *self.cbk_args[0], **self.cbk_kwargs[0])
                self._pos_unf.append(s.dump('pos'))

        # If trajectory is grandcanonical, we make sure all samples
        # have non-zero particles and raise an exception.
        for data in [self._pos]:
            if len(data) > 0:
                if 0 in [len(p) for p in self._pos]:
                    pass
                    #raise ValueError('cannot handle null samples in GC trajectory')

    def _setup_arrays_twobody(self):
        if len(self.cbk) <= 1:
            self._setup_arrays_onebody()
            self._pos_0 = self._pos
            self._pos_1 = self._pos
            return

        self._pos_0, self._pos_1 = [], []
        self._vel_0, self._vel_1 = [], []
        if 'pos' in self._phasespace or 'vel' in self._phasespace:
            for s in self.trajectory:
                s0 = self.cbk[0](s, *self.cbk_args[0], **self.cbk_kwargs[0])
                s1 = self.cbk[1](s, *self.cbk_args[1], **self.cbk_kwargs[1])
                if 'pos' in self._phasespace:
                    self._pos_0.append(s0.dump('pos'))
                    self._pos_1.append(s1.dump('pos'))
                if 'vel' in self._phasespace:
                    self._vel_0.append(s0.dump('vel'))
                    self._vel_1.append(s1.dump('vel'))

        # Dump unfolded positions if requested
        self._pos_unf_0, self._pos_unf_1 = [], []
        if 'pos-unf' in self._phasespace:
            for s in Unfolded(self.trajectory):
                s0 = self.cbk[0](s, *self.cbk_args[0], **self.cbk_kwargs[0])
                s1 = self.cbk[1](s, *self.cbk_args[1], **self.cbk_kwargs[1])
                self._pos_unf_0.append(s0.dump('pos'))
                self._pos_unf_1.append(s1.dump('pos'))

        # If trajectory is grandcanonical, we make sure all samples
        # have non-zero particles and raise an exception.
        for data in [self._pos_0, self._pos_1]:
            if 0 in [len(p) for p in data]:
                #raise ValueError('cannot handle null samples in GC trajectory')
                pass

    def compute(self):
        # Log
        if not self._need_update:
            log.info('skip %s (%s) for %s' % (self.short_name, self.tag, self.trajectory.filename))
            return
        
        log.debug('setup')
        from atooms.core.utils import Timer
        t = [Timer(), Timer()]
        t[0].start()
        self._setup_arrays()
        t[0].stop()
        log.debug('compute')
        t[1].start()
        self._compute()
        t[1].stop()
        log.info('computed %s %s for %s in %.1f sec [setup:%.0f%%, compute: %.0f%%]', self.description,
                 self.tag_description, self.trajectory.filename, t[0].wall_time + t[1].wall_time,
                 t[0].wall_time / (t[0].wall_time + t[1].wall_time) * 100,
                 t[1].wall_time / (t[0].wall_time + t[1].wall_time) * 100)
        try:
            self.analyze()
        except ImportError as e:
            log.warn('no analysis due to missing modules')
        return self.grid, self.value

    def _compute(self):
        pass

    def analyze(self):
        pass

    @property
    def _output_file(self):
        # Interpolate the output path string
        if self.output_path is None:
            filename = None
        else:
            # Add self. prefix to all attributes
            path = self.output_path.replace('{', '{0.')
            try:
                filename = path.format(self)
            except:
                print path
                raise
            # Strip unpleasant punctuation
            for punct in ['.', '_', '-']:
                filename = filename.replace(punct*2, punct)
                filename = filename.strip(punct)
        return filename

    def read(self):
        inp = open(self._output_file, 'r')
        x = numpy.loadtxt(inp, unpack=True)
        if len(x) == 3:
            self.grid, self.value = x
        elif len(x) == 2:
            self.grid, self.value = x
        else:
            self.grid, self.value = x[0:2]
            warnings.warn("Ignoring some columns in %s" % self._output_file)

        inp.close()

    def write(self):
        # TODO: it is probably the compute method that should be responsible for dumping grids appropriately
        if (isinstance(self.grid[0], list) or \
            isinstance(self.grid[0], numpy.ndarray)) and \
            len(self.grid) == 2:
            x = numpy.array(self.grid[0]).repeat(len(self.value[0]))
            y = numpy.array(self.grid[1] * len(self.grid[0]))
            z = numpy.array(self.value).flatten()
            dump = numpy.transpose(numpy.array([x, y, z]))
        else:
            dump = numpy.transpose(numpy.array([self.grid, self.value]))

        # Comment line
        if isinstance(self.variables, tuple) or isinstance(self.variables, list):
            columns = list(self.variables) + [self.description]
        else:
            columns = [self.variables] + [self.description]
        if len(self.tag_description) > 0:
            conj = 'of'
        else:
            conj = ''
        comments = _dump(title='%s %s %s' % (self.description, conj, self.tag_description),
                         columns=columns,
                         command='atooms-pp', version=__version__,
                         description=None, note=None,
                         parents=self.trajectory.filename,
                         inline=False)
        if not self.comments is None:
            comments += self.comments

        # Analyze results
        analysis = ""
        for x, f in self.results.iteritems():
            if f is not None:
                analysis += '# %s: %s\n' % (x, f)

        # Put it all together
        with open(self._output_file, 'w') as fh:
            fh.write(comments)
            if len(analysis) > 0:
                fh.write(analysis)
            numpy.savetxt(fh, dump, fmt="%g")
            fh.flush()

    def do(self):
        if not self._need_update:
            return
        self.compute()
        try:
            self.analyze()
        except ImportError as e:
            print 'Could not analyze due to missing modules, continuing...'
            print e.message
        self.write()

    def __call__(self):
        self.do()


