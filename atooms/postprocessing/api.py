"""Post processing API."""

import atooms.postprocessing as pp
from atooms.postprocessing.partial import Partial
from atooms.trajectory import Trajectory
from atooms.trajectory.decorators import change_species
from atooms.system.particle import distinct_species

from .helpers import linear_grid, logx_grid

_func_db = {'linear_grid': linear_grid, 'logx_grid': logx_grid, 'linear': linear_grid, 'logx': logx_grid}

def _get_trajectories(input_files, args):
    from atooms.trajectory import Sliced
    from atooms.core.utils import fractional_slice
    for input_file in input_files:
        with Trajectory(input_file, fmt=args['fmt']) as th:
            # Caching is useful for systems with multiple species but
            # it will increase the memory footprint. Use --no-cache to
            # disable it
            if not args['no_cache']:
                th.cache = True
            if args['species_layout'] is not None:
                th.register_callback(change_species, args['species_layout'])
            sl = fractional_slice(args['first'], args['last'], args['skip'], len(th))
            if th.block_size > 1:
                sl_start = (sl.start // th.block_size) * th.block_size if sl.start is not None else sl.start
                sl_stop = (sl.stop // th.block_size) * th.block_size if sl.stop is not None else sl.stop
                sl = slice(sl_start, sl_stop, sl.step)
            ts = Sliced(th, sl)
            yield ts

def _compat(args, fmt=None, species_layout=None):
    if 'first' not in args:
        args['first'] = None
    if 'last' not in args:
        args['last'] = None
    if 'skip' not in args:
        args['skip'] = None
    if 'fmt' not in args:
        args['fmt'] = fmt
    if 'species_layout' not in args:
        args['species_layout'] = species_layout
    if 'norigins' not in args:
        args['norigins'] = None
    if 'fast' not in args:
        args['fast'] = False
    if 'no_cache' not in args:
        args['no_cache'] = False
    if 'update' not in args:
        args['update'] = False
    if 'filter' not in args:
        args['filter'] = None
    if 'no_partial' not in args:
        args['no_partial'] = False

    # Implict option rules
    if args['filter'] is not None:
        args['no_partial'] = True

    return args

def gr(input_file, dr=0.04, grandcanonical=False, fmt=None,
       species_layout=None, *input_files, **global_args):
    """Radial distribution function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        th._grandcanonical = grandcanonical

        cf = pp.RadialDistributionFunction(th, dr=dr, norigins=global_args['norigins'])
        if global_args['filter'] is not None:
            cf = pp.Filter(cf, global_args['filter'])
        cf.do(update=global_args['update'])

        ids = distinct_species(th[0].particle)
        if len(ids) > 1 and not global_args['no_partial']:
            cf = Partial(pp.RadialDistributionFunction, ids, th, dr=dr, norigins=global_args['norigins'])
            cf.do(update=global_args['update'])

def sk(input_file, nk=20, dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30,
       species_layout=None, trajectory_field=None,
       field=None, *input_files, **global_args):
    """Structure factor"""
    global_args = _compat(global_args, species_layout=species_layout)
    if global_args['fast']:
        backend = pp.StructureFactorOpti
    else:
        backend = pp.StructureFactor

    for th in _get_trajectories([input_file] + list(input_files), global_args):
        backend(th, None, norigins=global_args['norigins'],
                trajectory_field=trajectory_field,
                field=field, kmin=kmin,
                kmax=kmax, nk=nk, dk=dk,
                ksamples=ksamples).do(update=global_args['update'])

        ids = distinct_species(th[0].particle)
        if len(ids) > 1 and trajectory_field is None and field is None:
            Partial(backend, ids, th, None,
                    norigins=global_args['norigins'],
                    kmin=kmin, kmax=kmax, nk=nk, dk=dk,
                    ksamples=ksamples).do(update=global_args['update'])

def ik(input_file, trajectory_radius=None, nk=20, dk=0.1, kmin=-1.0, kmax=15.0,
       ksamples=30, fmt=None, species_layout=None,
       *input_files, **global_args):
    """Spectral density"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if trajectory_radius is None:
            trajectory_radius = input_file
            pp.SpectralDensity(th, trajectory_radius,
                               kgrid=None, norigins=global_args['norigins'],
                               kmin=kmin, kmax=kmax, nk=nk, dk=dk,
                               ksamples=ksamples).do(update=global_args['update'])

def msd(input_file, time_target=-1.0, time_target_fraction=0.75,
        tsamples=30, sigma=1.0, func=linear_grid, grid=None, rmsd_target=-1.0,
        fmt=None, species_layout=None, *input_files, **global_args):
    """Mean square displacement"""
    if grid is not None and grid in _func_db:
        func = _func_db[grid]
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        dt = th.timestep
        if time_target > 0:
            t_grid = [0.0] + func(dt, min(th.total_time, time_target), tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(dt, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        #_cache, th.cache = th.cache, False
        ids = distinct_species(th[0].particle)
        #th.cache = _cache
        #print th.cache, id(th._cache)
        pp.MeanSquareDisplacement(th, tgrid=t_grid,
                                              norigins=global_args['norigins'],
                                              sigma=sigma, rmax=rmsd_target).do(update=global_args['update'])
        #print th.cache, th._cache
        if len(ids) > 1:
            #print th.cache, id(th._cache)
            Partial(pp.MeanSquareDisplacement, ids,
                    th, tgrid=t_grid, norigins=global_args['norigins'], sigma=sigma, rmax=rmsd_target).do(update=global_args['update'])

def vacf(input_file, time_target=-1.0, time_target_fraction=0.10,
         tsamples=30, func=linear_grid, fmt=None, species_layout=None,
         *input_files, **global_args):
    """Velocity autocorrelation function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, min(th.total_time, time_target), tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        pp.VelocityAutocorrelation(th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.VelocityAutocorrelation, ids, th,
                    t_grid, norigins=global_args['norigins']).do(update=global_args['update'])

def fkt(input_file, time_target=-1.0, time_target_fraction=0.75,
        tsamples=60, kmin=7.0, kmax=7.0, ksamples=1, dk=0.1, nk=100,
        tag_by_name=False, func=logx_grid, fmt=None,
        species_layout=None, *input_files, **global_args):
    """Total intermediate scattering function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        k_grid = linear_grid(kmin, kmax, ksamples)
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.IntermediateScattering, ids, th, k_grid, t_grid,
                    nk=nk, dk=dk).do(update=global_args['update'])

def fskt(input_file, time_target=-1.0, time_target_fraction=0.75,
         tsamples=60, kmin=7.0, kmax=8.0, ksamples=1, dk=0.1, nk=8,
         func=logx_grid, fmt=None, species_layout=None, total=False,
         *input_files, **global_args):
    """Self intermediate scattering function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        k_grid = linear_grid(kmin, kmax, ksamples)
        if total:
            pp.SelfIntermediateScattering(th, k_grid, t_grid,
                                                      nk, dk=dk, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.SelfIntermediateScattering, ids,
                    th, k_grid, t_grid, nk, dk=dk, norigins=global_args['norigins']).do(update=global_args['update'])

def chi4qs(input_file, tsamples=60, a=0.3, time_target=-1.0,
           time_target_fraction=0.75, fmt=None, species_layout=None,
           total=False, *input_files, **global_args):
    """Dynamic susceptibility of self overlap"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)

    if global_args['fast']:
        backend = pp.Chi4SelfOverlapOpti
    else:
        backend = pp.Chi4SelfOverlap

    for th in _get_trajectories([input_file] + list(input_files), global_args):
        func = logx_grid
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, min(th.total_time, time_target), tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        if total:
            backend(th, t_grid, a=a, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if not total and len(ids) > 1:
            Partial(backend, ids, th, t_grid, a=a, norigins=global_args['norigins']).do(update=global_args['update'])

def alpha2(input_file, time_target=-1.0, time_target_fraction=0.75,
           tsamples=60, func=logx_grid, fmt=None,
         species_layout=None, *input_files, **global_args):
    """Non-Gaussian parameter"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        pp.NonGaussianParameter(th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.NonGaussianParameter, ids, th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])

def qst(input_file, time_target=-1.0, time_target_fraction=0.75,
        tsamples=60, func=logx_grid, fmt=None,
        species_layout=None, *input_files, **global_args):
    """Self overlap correlation function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        pp.SelfOverlap(th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.SelfOverlap, ids, th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])

def qt(input_file, time_target=-1.0, time_target_fraction=0.75,
        tsamples=60, func=logx_grid, fmt=None,
        species_layout=None, *input_files, **global_args):
    """Collective overlap correlation function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        pp.CollectiveOverlap(th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(pp.CollectiveOverlap, ids, th, t_grid, norigins=global_args['norigins']).do(update=global_args['update'])

def ba(input_file, dtheta=4.0, grandcanonical=False, fmt=None,
       species_layout=None, *input_files, **global_args):
    """Bond-angle distribution"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        th._grandcanonical = grandcanonical

        cf = pp.BondAngleDistribution(th, dtheta=dtheta, norigins=global_args['norigins'])
        if global_args['filter'] is not None:
            cf = pp.Filter(cf, global_args['filter'])
        cf.do(update=global_args['update'])

        # ids = distinct_species(th[0].particle)
        # if len(ids) > 1 and not global_args['no_partial']:
        #     cf = Partial(pp.BondAngleDistribution, ids, th, dtheta=dtheta, norigins=global_args['norigins'])
        #     cf.do(update=global_args['update'])
