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

def filter_species(system, species):
    """Callback to filter particles by species.

    The input species can be an integer (particle id), a string
    (particle name), or None. In this latter case, all particles
    are returned.
    """
    s = copy.copy(system)
    if species is not None:
        if type(species) is not str:
            s.particle = [p for p in system.particle if p.id == species]
        else:
            s.particle = [p for p in system.particle if p.name == species]
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
