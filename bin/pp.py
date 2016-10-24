#!/usr/bin/env python

from atooms import postprocessing
from atooms.trajectory import Trajectory
from pyutils.utils import linear_grid, logx_grid

"""Test of a post processing runner."""

def msd(fname, msd_target=3.0, time_target=-1.0, t_samples=30, skip=1, sigma=1.0, func=linear_grid):
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
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    cf = postprocessing.VelocityAutocorrelation(trajectory, t_grid)
    cf.do()
    cf.do_species()
    #cf.do_dims()
    trajectory.close()

def fkt(fname, time_target=1e9, t_samples=60, k_min=7.0, k_max=7.0, k_samples=1, dk=0.1, tag_by_name=False, func=logx_grid):
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    k_grid = linear_grid(k_min, k_max, k_samples)
    cf = postprocessing.IntermediateScattering(trajectory, k_grid, t_grid)
    cf.do()
    cf.do_species(tag_by_name)
    trajectory.close()

def fskt(fname, time_target=1e9, t_samples=60, k_min=7.0, k_max=8.0, k_samples=1, dk=0.1, tag_by_name=False, func=logx_grid):
    trajectory = Trajectory(fname)
    t_grid = [0.0] + func(trajectory.timestep, time_target, t_samples)
    k_grid = linear_grid(k_min, k_max, k_samples)
    cf = postprocessing.SelfIntermediateScattering(trajectory, k_grid, t_grid)
    cf.do()
    cf.do_species(tag_by_name)
    trajectory.close()

def s4kt(fname, time_target=-1.0, t_samples=60, k_min=-1, k_max=4.0, k_samples=30, dk=0.1, tag_by_name=False, func=linear_grid):
    trajectory = Trajectory(fname)
    cf = postprocessing.SelfIntermediateScattering(trajectory, linear_grid(k_min, k_max, k_samples), [time_target])
    cf.do()
    cf.do_species(tag_by_name)
    trajectory.close()

# TODO: expose gr too

if __name__ == '__main__':

    import sys
    import time
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--tag-by-name', action='store_true', dest='tag_by_name',   default=False)
    parser.add_option('--gr',       action='store_true', dest='gr',       default=False)
    parser.add_option('--msd',      action='store_true', dest='msd',      default=False)
    parser.add_option('--msd-samples', action='store',   dest='msd_samples', default=30,   type='int')
    parser.add_option('--msd-sigma', action='store',      dest='msd_sigma', default=1.0, type='float')
    parser.add_option('--msd-tmax', action='store',      dest='msd_tmax', default=-1.0, type='float')
    parser.add_option('--msd-target',action='store',     dest='msd_target',default=3.0, type='float')
    parser.add_option('--vacf',     action='store_true', dest='vacf',     default=False)
    parser.add_option('--vacf-samples', action='store',   dest='vacf_samples', default=60,   type='int')
    parser.add_option('--vacf-tmax', action='store',      dest='vacf_tmax', default=1.0, type='float')
    parser.add_option('--fkt',         action='store_true', dest='fkt',     default=False)
    parser.add_option('--fkt-tsamples', action='store',      dest='fkt_samples', default=60,   type='int')
    parser.add_option('--fkt-tmax',     action='store',      dest='fkt_tmax', default=1e9, type='float')
    parser.add_option('--fkt-kmin',     action='store',      dest='fkt_kmin', default=7.0, type='float')
    parser.add_option('--fkt-kmax',     action='store',      dest='fkt_kmax', default=8.0, type='float')
    parser.add_option('--fkt-dk',       action='store',      dest='fkt_dk',   default=0.1, type='float')
    parser.add_option('--fkt-ksamples', action='store',      dest='fkt_ksamples', default=1, type='int')
    parser.add_option('--fkt-part',     action='store',      dest='fkt_part', default='both', type='str')
    (opts, args) = parser.parse_args()

    def gimmefile(what):
        sys.stdout.write("%s... " % what)
        sys.stdout.flush()
        return fname + '.pp.%s'

    for fname in args:
        t1 = time.time()
        sys.stdout.write("# [pp.py] analysing %s : getting " % fname)

        if opts.msd:
            sys.stdout.write("msd... ")
            sys.stdout.flush()
            msd(fname, msd_target=opts.msd_target, t_samples=opts.msd_samples, time_target=opts.msd_tmax)
        if opts.fkt:
            sys.stdout.write("fkt... ")
            sys.stdout.flush()
            if opts.fkt_part == 'both' or opts.fkt_part == 'self':
                fskt(fname, opts.fkt_tmax, opts.fkt_samples, opts.fkt_kmin, opts.fkt_kmax, opts.fkt_ksamples, opts.fkt_dk, opts.tag_by_name) 
            if opts.fkt_part == 'both' or opts.fkt_part == 'total':
                fkt(fname, opts.fkt_tmax, opts.fkt_samples, opts.fkt_kmin, opts.fkt_kmax, opts.fkt_ksamples, opts.fkt_dk, opts.tag_by_name) 
        if opts.vacf:
            sys.stdout.write("vacf... ")
            sys.stdout.flush()
            vacf(fname, opts.vacf_tmax, opts.vacf_samples)
        t2 = time.time()
        sys.stdout.write("# done in %.1f s\n" % (t2-t1))
        
