"""Post processing API."""

import postprocessing
from atooms.trajectory import Trajectory
from pyutils.utils import linear_grid, logx_grid


def msd(fname, msd_target=3.0, time_target=-1.0, t_samples=30, skip=1, sigma=1.0, func=linear_grid):
    """Mean square displacement."""
    trajectory = Trajectory(fname)
    dt = trajectory.timestep
    if time_target < 0.0:
        t_grid = [0.0] + func(dt, trajectory.time_when_msd_is(msd_target), t_samples)
    else:
        t_grid = [0.0] + func(dt, min(trajectory.steps[-1]*dt, time_target), t_samples)
    cf = postprocessing.MeanSquareDisplacement(trajectory, t_grid, skip, sigma)
    #TODO: optimize avoid doing the species calculation is nsp=1
    #TODO: switch to enable xyz
    cf.do()
    cf.do_species()
    cf.do_dims()
    trajectory.close()

def vacf(fname, time_target=1.0, t_samples=30, func=linear_grid):
    """Velocity autocorrelation function."""
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    cf = postprocessing.VelocityAutocorrelation(trajectory, t_grid)
    cf.do()
    cf.do_species()
    #cf.do_dims()
    trajectory.close()

def fkt(fname, time_target=1e9, t_samples=60, k_min=7.0, k_max=7.0, k_samples=1, dk=0.1, tag_by_name=False, func=logx_grid):
    """Total intermediate scattering function."""
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    k_grid = linear_grid(k_min, k_max, k_samples)
    cf = postprocessing.IntermediateScattering(trajectory, k_grid, t_grid)
    cf.do()
    cf.do_species(tag_by_name)
    trajectory.close()

def fskt(fname, time_target=1e9, t_samples=60, k_min=7.0, k_max=8.0, k_samples=1, dk=0.1, tag_by_name=False, func=logx_grid):
    """Self intermediate scattering function."""
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    k_grid = linear_grid(k_min, k_max, k_samples)
    cf = postprocessing.SelfIntermediateScattering(trajectory, k_grid, t_grid)
    cf.do()
    cf.do_species(tag_by_name)
    trajectory.close()
