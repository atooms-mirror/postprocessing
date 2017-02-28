"""Post processing API."""

import postprocessing
from postprocessing.partial import Partial
from atooms.trajectory import Trajectory
from atooms.system.particle import species
from .helpers import linear_grid, logx_grid

# TODO; can we set trajectory format globally?

def gr(input_file, grandcanonical=False, fmt=None):
    """Radial distribution function."""
    with Trajectory(input_file, fmt=fmt) as th:
        th._grandcanonical = grandcanonical
        postprocessing.RadialDistributionFunction(th).do()
        # We only do partial correlation functions in mixtures
        ids = species(th[0].particle)
        if len(ids) > 1:
            Partial(postprocessing.RadialDistributionFunction, [1, 2], th).do()

def sk(input_file, nk=20, dk=0.1, kmin=-1.0, kmax=15.0, ksamples=30, norigins=-1,
       grandcanonical=False, fmt=None):
    """Structure factor."""    
    with Trajectory(input_file, fmt=fmt) as th:
        k_grid = linear_grid(kmin, kmax, ksamples)
        postprocessing.StructureFactor(th, k_grid, norigins=norigins).do()
        # We only do partial correlation functions in mixtures
        ids = species(th[0].particle)
        if len(ids) > 1:
            Partial(postprocessing.StructureFactor, ids, th, k_grid).do()

def msd(input_file, msd_target=3.0, time_target=-1.0, t_samples=30,
        norigins=50, sigma=1.0, func=linear_grid, fmt=None):
    """Mean square displacement."""
    with Trajectory(input_file, fmt=fmt) as th:
        dt = th.timestep
        # TODO: restore this and get rid of relic in trjectory class
        # if time_target < 0.0:
        #     t_grid = [0.0] + func(dt, th.time_when_msd_is(msd_target), t_samples)
        # else:
        #t_grid = [0.0] + func(dt, min(th.steps[-1]*dt, time_target), t_samples)
        t_grid = [0.0] + func(dt, th.steps[-1]*dt, t_samples)
        # TODO: onebody correlations should be filtered by single species, not pair
        for species in [1, 2]:
            cf = Partial(postprocessing.MeanSquareDisplacement, [species],
                         th, tgrid=t_grid, norigins=norigins, sigma=sigma)
            cf.do()
        cf = postprocessing.MeanSquareDisplacement(th, tgrid=t_grid, norigins=norigins, sigma=sigma)
        cf.do()

def vacf(input_file, time_target=1.0, t_samples=30, func=linear_grid, fmt=None):
    """Velocity autocorrelation function."""
    with Trajectory(input_file, fmt=fmt) as trajectory:
        t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
        for species in [1, 2]:
            cf = Partial(postprocessing.VelocityAutocorrelation, [species], trajectory, t_grid)
            cf.do()
        cf = postprocessing.VelocityAutocorrelation(trajectory, t_grid)
        cf.do()

def fkt(input_file, time_target=1e9, tsamples=60, kmin=7.0, kmax=7.0,
        ksamples=1, dk=0.1, tag_by_name=False, func=logx_grid, fmt=None):
    """Total intermediate scattering function."""
    with Trajectory(input_file, fmt=fmt) as trajectory:
        t_grid = [0.0] + func(trajectory.timestep, time_target, tsamples)
        k_grid = linear_grid(kmin, kmax, ksamples)
        cf = postprocessing.IntermediateScattering(trajectory, k_grid, t_grid)
        cf.do()
        cf.do_species(tag_by_name)

def fskt(input_file, time_target=1e9, tsamples=60, kmin=7.0, kmax=8.0,
         ksamples=1, dk=0.1, tag_by_name=False, func=logx_grid, fmt=None):
    """Self intermediate scattering function."""
    with Trajectory(input_file, fmt=fmt) as trajectory:
        t_grid = [0.0] + func(trajectory.timestep, time_target, tsamples)
        k_grid = linear_grid(kmin, kmax, ksamples)
        cf = postprocessing.SelfIntermediateScattering(trajectory, k_grid, t_grid)
        cf.do()
        cf.do_species(tag_by_name)
