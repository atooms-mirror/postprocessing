"""Post processing API."""

import postprocessing
from postprocessing.partial import Partial
from atooms.core.utils import setup_logging
from atooms.trajectory import Trajectory
from atooms.trajectory.decorators import filter_species
from atooms.trajectory.utils import time_when_msd_is
from atooms.system.particle import distinct_species
from .helpers import linear_grid, logx_grid

setup_logging('postprocessing', level=20)
setup_logging('atooms', level=20)

def gr(input_file, grandcanonical=False, fmt=None, show=False, norigins=-1):
    """Radial distribution function."""
    with Trajectory(input_file, fmt=fmt) as th:
        th._grandcanonical = grandcanonical
        postprocessing.RadialDistributionFunction(th, norigins=norigins).do(show)
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.RadialDistributionFunction, ids, th).do(show)

def sk(input_file, nk=20, dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30,
       norigins=-1, species=None, grandcanonical=False, fmt=None,
       trajectory_field=None, field=None):
    """Structure factor."""
    with Trajectory(input_file, fmt=fmt) as th:
        if species is not None:
            th.register_callback(filter_species, int(species))
        ids = distinct_species(th[-1].particle)
        postprocessing.StructureFactor(th, None, norigins=norigins,
                                       trajectory_field=trajectory_field,
                                       field=field, kmin=kmin,
                                       kmax=kmax, nk=nk,
                                       ksamples=ksamples).do()
        if len(ids) > 1 and trajectory_field is None:
            Partial(postprocessing.StructureFactor, ids, th, None,
                    norigins=norigins, kmin=kmin,
                    kmax=kmax, nk=nk,
                    ksamples=ksamples).do()

def ik(input_file, trajectory_radius=None, nk=20, dk=0.1, kmin=-1.0,
       kmax=15.0, ksamples=30, norigins=-1, species=None,
       grandcanonical=False, fmt=None):
    """Spectral density,"""
    if trajectory_radius is None:
        trajectory_radius = input_file
    with Trajectory(input_file, fmt=fmt) as th:
        if species is not None:
            th.register_callback(filter_species, int(species))
        ids = distinct_species(th[-1].particle)
        postprocessing.SpectralDensity(th, trajectory_radius,
                                       kgrid=None, norigins=norigins,
                                       kmin=kmin, kmax=kmax, nk=nk,
                                       ksamples=ksamples).do()

def msd(input_file, rmsd_target=-1.0, time_target=-1.0,
        time_target_fraction=-1.0, tsamples=30, norigins=50, sigma=1.0,
        func=linear_grid, fmt=None):
    """Mean square displacement."""
    with Trajectory(input_file, fmt=fmt) as th:
        dt = th.timestep
        if rmsd_target > 0:
            t_grid = [0.0] + func(dt, time_when_msd_is(th, rmsd_target**2),
                                  tsamples)
        else:
            if time_target > 0:
                t_grid = [0.0] + func(dt, min(th.total_time,
                                              time_target), tsamples)
            elif time_target_fraction > 0:
                t_grid = [0.0] + func(dt, time_target_fraction*th.total_time,
                                      tsamples)
            else:
                t_grid = [0.0] + func(dt, th.steps[-1]*dt, tsamples)
        ids = distinct_species(th[-1].particle)
        postprocessing.MeanSquareDisplacement(th, tgrid=t_grid, norigins=norigins, sigma=sigma).do()
        if len(ids) > 1:
            Partial(postprocessing.MeanSquareDisplacement, ids,
                    th, tgrid=t_grid, norigins=norigins, sigma=sigma).do()

def vacf(input_file, time_target=1.0, tsamples=30, func=linear_grid, fmt=None):
    """Velocity autocorrelation function."""
    with Trajectory(input_file, fmt=fmt) as th:
        t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        postprocessing.VelocityAutocorrelation(th, t_grid).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.VelocityAutocorrelation, ids, th, t_grid).do()

def fkt(input_file, time_target=1e9, tsamples=60, kmin=7.0, kmax=7.0,
        ksamples=1, dk=0.1, tag_by_name=False, func=logx_grid, fmt=None):
    """Total intermediate scattering function."""
    with Trajectory(input_file, fmt=fmt) as th:
        t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        k_grid = linear_grid(kmin, kmax, ksamples)
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(postprocessing.IntermediateScattering, ids, th, k_grid, t_grid).do()

def fskt(input_file, time_target=1e9, tsamples=60, kmin=7.0, kmax=8.0,
         ksamples=1, dk=0.1, nk=8, skip=1, tag_by_name=False, func=None, fmt=None):
    """Self intermediate scattering function."""
    with Trajectory(input_file, fmt=fmt) as th:
        if func is None:
            func = logx_grid
            t_grid = [0.0] + func(th.timestep, min(th.times[-1], time_target), tsamples)
        else:
            t_grid = [th.timestep*i for i in th.steps]
        k_grid = linear_grid(kmin, kmax, ksamples)
        postprocessing.SelfIntermediateScattering(th, k_grid, t_grid, nk, dk=dk, skip=skip).do()
        ids = distinct_species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.SelfIntermediateScattering, ids, th, k_grid, t_grid, nk, dk=dk, skip=skip).do()

def chi4qs(input_file, tsamples=60, a=0.3, fmt=None):
    """Dynamic susceptibility of self overlap."""
    with Trajectory(input_file, fmt=fmt) as th:
        func = logx_grid
        time_target = th.total_time * 0.75
        t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        ids = distinct_species(th[0].particle)
        if len(ids) > 1:
            Partial(postprocessing.Chi4SelfOverlap, ids, th, t_grid, a=a).do()
        else:
            postprocessing.Chi4SelfOverlap(th, t_grid, a=a).do()
