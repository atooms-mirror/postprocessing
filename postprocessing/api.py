"""Post processing API."""

import postprocessing
from postprocessing.partial import Partial
from atooms.utils import setup_logging
from atooms.trajectory import Trajectory
from atooms.trajectory.decorators import filter_id
from atooms.trajectory.utils import time_when_msd_is
from atooms.system.particle import species
from .helpers import linear_grid, logx_grid

setup_logging('postprocessing', level=20)
setup_logging('atooms', level=20)

def gr(input_file, grandcanonical=False, fmt=None, show=False):
    """Radial distribution function."""
    with Trajectory(input_file, fmt=fmt) as th:
        th._grandcanonical = grandcanonical
        postprocessing.RadialDistributionFunction(th).do(show)
        ids = species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.RadialDistributionFunction, ids, th).do(show)

def sk(input_file, nk=20, dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30, norigins=-1, include_id=None,
       grandcanonical=False, fmt=None):
    """Structure factor."""    
    with Trajectory(input_file, fmt=fmt) as th:
        k_grid = linear_grid(kmin, kmax, ksamples)
        if include_id is not None:
            th.register_callback(filter_id, int(include_id))
        ids = species(th[-1].particle)
        postprocessing.StructureFactor(th, k_grid, norigins=norigins).do()
        if len(ids) > 1:
            Partial(postprocessing.StructureFactor, ids, th, k_grid).do()

def msd(input_file, rmsd_target=-1.0, time_target=-1.0, t_samples=30,
        norigins=50, sigma=1.0, func=linear_grid, fmt=None):
    """Mean square displacement."""
    with Trajectory(input_file, fmt=fmt) as th:
        dt = th.timestep
        if rmsd_target > 0:
            t_grid = [0.0] + func(dt, time_when_msd_is(th, rmsd_target**2), t_samples)
        else:
            if time_target > 0:
                t_grid = [0.0] + func(dt, min(th.steps[-1]*dt, time_target), t_samples)
            else:
                t_grid = [0.0] + func(dt, th.steps[-1]*dt, t_samples)
        ids = species(th[-1].particle)
        postprocessing.MeanSquareDisplacement(th, tgrid=t_grid, norigins=norigins, sigma=sigma).do()
        if len(ids) > 1:
            Partial(postprocessing.MeanSquareDisplacement, ids,
                    th, tgrid=t_grid, norigins=norigins, sigma=sigma).do()

def vacf(input_file, time_target=1.0, t_samples=30, func=linear_grid, fmt=None):
    """Velocity autocorrelation function."""
    with Trajectory(input_file, fmt=fmt) as th:
        t_grid = [0.0] + func(th.timestep, time_target, t_samples)
        postprocessing.VelocityAutocorrelation(th, t_grid).do()
        ids = species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.VelocityAutocorrelation, ids, th, t_grid).do()

def fkt(input_file, time_target=1e9, tsamples=60, kmin=7.0, kmax=7.0,
        ksamples=1, dk=0.1, tag_by_name=False, func=logx_grid, fmt=None):
    """Total intermediate scattering function."""
    with Trajectory(input_file, fmt=fmt) as th:
        t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        k_grid = linear_grid(kmin, kmax, ksamples)
        ids = species(th[0].particle)
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
        ids = species(th[-1].particle)
        if len(ids) > 1:
            Partial(postprocessing.SelfIntermediateScattering, ids, th, k_grid, t_grid, nk, dk=dk, skip=skip).do()

def chi4qs(input_file, tsamples=60, a=0.3, fmt=None):
    """Dynamic susceptibility of self overlap."""
    with Trajectory(input_file, fmt=fmt) as th:
        func = logx_grid
        time_target = th.time_total * 0.75
        t_grid = [0.0] + func(th.timestep, time_target, tsamples)
        ids = species(th[0].particle)
        if len(ids) > 1:
            Partial(postprocessing.Chi4SelfOverlap, ids, th, t_grid, a).do()
        else:
            postprocessing.Chi4SelfOverlap(th, t_grid, a).do()
            
