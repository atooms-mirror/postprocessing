"""Post processing API."""

import atooms.postprocessing as postprocessing
from atooms.postprocessing.partial import Partial
from atooms.trajectory import Trajectory
from atooms.trajectory.decorators import change_species
from atooms.system.particle import distinct_species

from .helpers import linear_grid, logx_grid


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

def _compat(args, fmt, species_layout=None):
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
    return args

def gr(input_file, dr=0.04, grandcanonical=False, fmt=None, species_layout=None,
       *input_files, **global_args):
    """Radial distribution function"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        th._grandcanonical = grandcanonical
        postprocessing.RadialDistributionFunction(th, dr=dr, norigins=global_args['norigins']).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.RadialDistributionFunction, ids, th, dr=dr, norigins=global_args['norigins']).do()

def sk(input_file, nk=20, dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30,
       species_layout=None, fmt=None, trajectory_field=None,
       field=None, *input_files, **global_args):
    """Structure factor"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    if global_args['fast']:
        backend = postprocessing.StructureFactorOpti
    else:
        backend = postprocessing.StructureFactor

    for th in _get_trajectories([input_file] + list(input_files), global_args):
        ids = distinct_species(th[-1].particle)
        backend(th, None, norigins=global_args['norigins'],
                trajectory_field=trajectory_field,
                field=field, kmin=kmin,
                kmax=kmax, nk=nk, dk=dk,
                ksamples=ksamples).do()
        if len(ids) > 1 and trajectory_field is None:
            Partial(backend, ids, th, None,
                    norigins=global_args['norigins'], kmin=kmin,
                    kmax=kmax, nk=nk, dk=dk,
                    ksamples=ksamples).do()

def ik(input_file, trajectory_radius=None, nk=20, dk=0.1, kmin=-1.0, kmax=15.0,
       ksamples=30, fmt=None, species_layout=None,
       *input_files, **global_args):
    """Spectral density"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        if trajectory_radius is None:
            trajectory_radius = input_file
            postprocessing.SpectralDensity(th, trajectory_radius,
                                           kgrid=None, norigins=global_args['norigins'],
                                           kmin=kmin, kmax=kmax, nk=nk, dk=dk,
                                           ksamples=ksamples).do()

def msd(input_file, time_target=-1.0, time_target_fraction=0.75,
        tsamples=30, sigma=1.0, func=linear_grid, rmsd_target=-1.0,
        fmt=None, species_layout=None, *input_files, **global_args):
    """Mean square displacement"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)
    for th in _get_trajectories([input_file] + list(input_files), global_args):
        dt = th.timestep
        if time_target > 0:
            t_grid = [0.0] + func(dt, min(th.total_time, time_target), tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(dt, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        ids = distinct_species(th[-1].particle)
        postprocessing.MeanSquareDisplacement(th, tgrid=t_grid,
                                              norigins=global_args['norigins'],
                                              sigma=sigma, rmax=rmsd_target).do()
        if len(ids) > 1:
            Partial(postprocessing.MeanSquareDisplacement, ids,
                    th, tgrid=t_grid, norigins=global_args['norigins'], sigma=sigma, rmax=rmsd_target).do()

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
        postprocessing.VelocityAutocorrelation(th, t_grid,
                                               norigins=global_args['norigins']).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.VelocityAutocorrelation, ids, th,
                    t_grid, norigins=global_args['norigins']).do()

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
            Partial(postprocessing.IntermediateScattering, ids, th, k_grid, t_grid,
                    nk=nk, dk=dk).do()

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
            postprocessing.SelfIntermediateScattering(th, k_grid, t_grid,
                                                      nk, dk=dk, norigins=global_args['norigins']).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.SelfIntermediateScattering, ids,
                    th, k_grid, t_grid, nk, dk=dk, norigins=global_args['norigins']).do()

def chi4qs(input_file, tsamples=60, a=0.3, time_target=-1.0,
           time_target_fraction=0.75, fmt=None, species_layout=None,
           total=False, *input_files, **global_args):
    """Dynamic susceptibility of self overlap"""
    global_args = _compat(global_args, fmt=fmt, species_layout=species_layout)

    if global_args['fast']:
        backend = postprocessing.Chi4SelfOverlapOpti
    else:
        backend = postprocessing.Chi4SelfOverlap

    for th in _get_trajectories([input_file] + list(input_files), global_args):
        func = logx_grid
        if time_target > 0:
            t_grid = [0.0] + func(th.timestep, min(th.total_time, time_target), tsamples)
        elif time_target_fraction > 0:
            t_grid = [0.0] + func(th.timestep, time_target_fraction*th.total_time, tsamples)
        else:
            t_grid = None
        if total:
            backend(th, t_grid, a=a, norigins=global_args['norigins']).do()
        ids = distinct_species(th[0].particle)
        if not total and len(ids) > 1:
            Partial(backend, ids, th, t_grid, a=a, norigins=global_args['norigins']).do()

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
        postprocessing.NonGaussianParameter(th, t_grid, norigins=global_args['norigins']).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.NonGaussianParameter, ids, th, t_grid, norigins=global_args['norigins']).do()
