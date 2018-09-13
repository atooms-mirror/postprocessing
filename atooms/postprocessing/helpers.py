import copy


def linear_grid(min,max,delta):
    """Linear grid."""
    if type(delta) is int:
        n = delta
        if n > 1:
            delta = (max - min) / (n-1)
        else:
            delta = 0.0
    else:
        n = int((max-min)/delta)+1
    list = [min+i*delta for i in range(n)]
    return list


def logx_grid(x1, x2, n):
    """Create a list of n numbers in logx scale from x1 to x2."""
    # the shape if a*x^n. if n=0 => a=x1, if n=N => x1*x^N=x2
    if x1 > 0:
        xx = (x2/x1)**(1.0/n)
        return [x1] + [x1 * xx**(i+1) for i in range(1,n)]
    else:
        xx = (x2)**(1.0/n)
        return [x1] + [xx**(i+1)-1 for i in range(1,n)]


def ifabsmm(x, f):
    """Interpolated absolute maximum."""

    def _vertex_parabola(a,b,c):
        """Returns the vertex (x,y) of a parabola of the type a*x**2 + b*x + c."""
        return -b/(2*a), - (b**2 - 4*a*c) / (4*a)

    def _parabola_3points(x1,y1,x2,y2,x3,y3):
        """Parabola through 3 points."""
        delta = (x1 - x2)*(x1 - x3)*(x2 - x3)
        a     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / delta
        b     = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / delta
        c     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / delta
        return a, b, c

    # First uninterpolated minima and maxima
    imin, imax = f.index(min(f)), f.index(max(f))
    # Then perform parabolic interpolation
    ii = []
    for i in [imin, imax]:
        i1 = i-1
        i2 = i
        i3 = i+1
        a, b, c = _parabola_3points(x[i1], f[i1], x[i2], f[i2], x[i3], f[i3])
        ii.append(_vertex_parabola(a, b, c))
    return ii[0], ii[1]


def linear_fit(xdata, ydata):
    """
    Linear regression.

    Expressions as in Wikipedia (https://en.wikipedia.org/wiki/Simple_linear_regression)
    """
    import numpy
    from math import sqrt
    n = len(ydata)
    dof = n - 2
    sx = numpy.sum(xdata)
    sy = numpy.sum(ydata)
    sxy = sum(xdata * ydata)
    sxx = sum(xdata**2)
    syy = sum(ydata**2)

    a = (n * sxy - sx * sy) / (n * sxx - sx**2)
    b = sy / n - a * sx / n
    s = (n*syy - sy**2 - a**2 * (n*sxx - sx**2)) / (n*dof)
    sa = n * s / (n*sxx - sx**2)
    sb = sa * sxx / n

    return a, b, sqrt(sa), sqrt(sb)


def feqc(x, f, fstar):
    """
    Find first root of f=f(x) for data sets.

    Given two lists x and f, it returns the value of xstar for which
    f(xstar) = fstar. Raises an ValueError if no root is found.
    """
    s = f[0] - fstar
    for i in range(min(len(x), len(f))):
        if (f[i] - fstar) * s < 0.0:
            # Linear interpolation
            dxf   = (f[i] - f[i-1]) / (x[i] - x[i-1])
            xstar = x[i-1] + (fstar - f[i-1]) / dxf
            istar = i
            return xstar, istar

    # We get to the end and cannot find the root
    return None, None


def filter_species(system, species):
    """Callback to filter particles by species.

    The input species can be an integer (particle id), a string
    (particle name), or None. In this latter case, all particles
    are returned.
    """
    s = copy.copy(system)
    if species is not None:
        s.particle = [p for p in system.particle if p.species == species]
    return s

def filter_all(system):
    s = copy.copy(system)
    s.particle = [p for p in system.particle]
    return s

def adjust_skip(trajectory, n_origin=-1):
    """ Utility function to set skip so as to keep computation time under control """
    # TODO: We should also adjust it for Npart
    if trajectory.block_size > 1:
        return trajectory.block_size
    else:
        if n_origin > 0:
            return max(1, int(len(trajectory.steps) / float(n_origin)))
        else:
            return 1

def setup_t_grid(trajectory, t_grid):
    def templated(entry, template, keep_multiple=False):
        """Filter a list of entries so as to best match an input
        template. Lazy, slow version O(N*M). Ex.:
        entry=[1,2,3,4,5,10,20,100], template=[1,7,12,80] should
        return [1,5,10,100].
        """
        match = [min(entry, key=lambda x: abs(x-t)) for t in template]
        if not keep_multiple:
            match = list(set(match))
        return sorted(match)

    # First get all possible time differences
    steps = trajectory.steps
    off_samp = {}
    for off in range(trajectory.block_size):
        for i in range(off, len(steps)-off):
            if not steps[i] - steps[off] in off_samp:
                off_samp[steps[i] - steps[off]] = (off, i-off)

    # Retain only those pairs of offsets and sample
    # difference that match the desired input. This is the grid
    # used internally to calculate the time correlation function.
    i_grid = set([int(round(t/trajectory.timestep)) for t in t_grid])
    offsets = [off_samp[t] for t in templated(sorted(off_samp.keys()), sorted(i_grid))]
    # TODO: add this as a test
    # check = []
    # for off, i in offsets:
    #     for i0 in xrange(off, len(trajectory)-i, trajectory.block_size):
    #         check.append(trajectory.steps[i0+i] - trajectory.steps[i0])
    # print sorted(set(check)), sorted(dt)
    return offsets


# Add dump() from medepy for portability.
# This should be dropped in the future.
def _dump(title, columns=None, command=None, version=None,
          description=None, note=None, parents=None, inline=False,
          comment='# ', extra_fields=None):
    """
    Return a string of comments filled with metadata.
    """
    import datetime
    import os
    date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if columns is not None:
        columns_string  = ', '.join(columns)

    # Checksums of parent files
    if parents is not None:
        try:
            import md5
            # Make sure parents is list
            if not hasattr(parents, '__iter__'):
                parents = [parents]
            # Compute checksum
            checksums = []
            size_limit = 1e9
            if max([os.path.getsize(f) for f in parents]) < size_limit:
                for parentpath in parents:
                    tag = md5.md5(open(parentpath).read()).hexdigest()
                    checksums.append(tag)
                checksums = ', '.join(checksums)
            else:
                checksums = None
            # Convert to string
            parents = ', '.join([os.path.basename(p) for p in parents])
        except ImportError:
            checksums = None

    metadata = [('title', title),
                ('columns', columns_string),
                ('date', date),
                ('command', command),
                ('version', version),
                ('parents', parents),
                ('checksums', checksums),
                ('description', description),
                ('note', note)]

    if extra_fields is not None:
        metadata += extra_fields

    if inline:
        fmt = '{}: {};'
        txt = comment + ' '.join([fmt.format(key, value) for key,
                                  value in metadata if value is not None])
    else:
        txt = ''
        for key, value in metadata:
            if value is not None:
                txt += comment + '{}: {}\n'.format(key, value)
    return txt
