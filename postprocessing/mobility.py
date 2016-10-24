#!/usr/bin/env python
# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich


import os, sys
from math import acos, sqrt
import numpy
from atooms import trajectory as trj
import .realspace

def mobility(f, fn, tau):

    def kernel(t0, t1):
        return realspace.square_displacement(t1.dump('pos'), t0.dump('pos'))

    #tau=float(sys.argv[2])
    to = trj.Unfolded(trj.Trajectory(f))
    t = to.convert(trj.TrajectoryHDF5, tag='unf')
    #tm = trj.TrajectoryXYZ(f + '.mobility.tau%g' % tau, 'w')
    tm = trj.TrajectoryXYZ(f + '.mobility', 'w')
    tm.fmt.remove('velocity')
    # TODO: tm.metadata = t.metadata or even more compact
    tm.write_initial_state(t.read(0))

    if t.block_period == 1:
        i0 = int(tau / (t.timestep * (t.steps[1]-t.steps[0])))
        i0 = min(max(1, i0), len(t.steps)/2)
    else:
        raise ValueError('cannot do with block period > 1')
        
    # Trajectory from which we want to get the reference samples (could be ignored) typically structure calculation
    tn = trj.TrajectoryXYZ(fn)
    # Find all samples that match the reference steps
    iall = [t.steps.index(itarget) for itarget in tn.steps]
    for i in iall:
        if i+i0 >= len(t.samples):
            break
        dr = kernel(t.read(i), t.read(i+i0))
        s = t.read(i)
        for pi, dri in zip(s.particle, dr):
            pi.tag = dri
        tm.write_sample(s, t.steps[i], i)

        # for j, dri in enumerate(dr):
        #     if dri>10.0:
        #         for ii in range(i0):
        #             r = t.read(ii+i).particle[j].position
        #             print r[0], r[1], r[2]
        #         tm.close()
        #         t.close()
        #         os.system('rm -f ' + t.filename)
        #         return
        
    tm.close()
    t.close()
    os.system('rm -f ' + t.filename)

def angle(f, dr_max):

    def kernel(s0, s1):
        return s1.dump('pos') - s0.dump('pos')

    def angle(u, v):
        unorm = sum(u**2)
        vnorm = sum(v**2)
        return acos(numpy.dot(u, v) / sqrt(unorm*vnorm))

    t = trj.Unfolded(trj.Trajectory(f))

    # Find all samples that match the reference steps
    #iall = [t.steps.index(itarget) for itarget in tn.steps]
    npart = len(t.read(0).particle)
    from copy import deepcopy

    def jumps():
        """Compute angles at jumps"""
        dr = [[] for i in range(npart)]
        s0 = deepcopy(t.read(0).particle)
        ilast = [0] * npart
        for i0 in range(1000): #t.samples:
            #print i0, len(t.samples)
            s1 = deepcopy(t.read(i0).particle)
            for j in range(npart):
                p0 = deepcopy(s0[j].position)
                p1 = deepcopy(s1[j].position)
                if sqrt(sum((p1-p0)**2)) > dr_max:
                    dr[j].append(p1-p0)
                    s0[j] = deepcopy(s1[j])
                    #print 'appended', j, i0-ilast[j], sqrt(sum((p1-p0)**2)), p1-p0
                    ilast[j] = i0

        for dr_part in dr:
            for dr_0, dr_1 in zip(dr_part[:-1], dr_part[1:]):
                theta = angle(dr_0, dr_1)
                print theta / 3.1415 * 180.0 #dr_0, dr_1, theta

    def plain():
        """Compute angles at fixed times i0+i and i0+2*i"""
        if t.block_period == 1:
            i = int(tau / (t.timestep * (t.steps[1]-t.steps[0])))
            i = min(max(1, i), len(t.steps)/2)
        else:
            raise ValueError('cannot do with block period > 1')
        for i0 in t.samples[0::2*i]:
            if i0 + 2*i >= len(t.samples):
                break
            s0 = deepcopy(t.read(i0))
            s1 = deepcopy(t.read(i0+1*i))
            s2 = deepcopy(t.read(i0+2*i))

            for j in range(npart):
                p0 = s0.particle[j].position
                p1 = s1.particle[j].position
                p2 = s2.particle[j].position
                theta = angle(p1-p0, p2-p1)
                print theta, sqrt(sum((p2-p0)**2))

    jumps()

    t.close()

def displacement(ti, tau, out_dir, movie=False):
    from copy import deepcopy
    from atooms.plugins import visualize
    def kernel(s0, s1):
        return s1.dump('pos') - s0.dump('pos')

    # To avoid unfolding issues, we open two copies of the trajectory
    # TODO: this does not work on a decorated trajectory, it passes through new() not init()
    ti0 = ti.__class__(ti.filename, 'r')
    t = trj.Unfolded(ti)
    t0 = trj.Unfolded(ti0)
    if tau > 0:
        i = tau
    else:
        i = t.block_period

    # if t.block_period == 1:
    #     i = int(tau / (t.timestep * (t.steps[1]-t.steps[0])))
    #     i = min(max(1, i), len(t.steps)/2)
    # else:
    #     raise ValueError('cannot do with block period > 1')

    k = 1
    ijump = t.block_period

    print '# Time interval     :', 0.2* t.timestep*(t.steps[i]-t.steps[0])
    print '# Time origin delta :', 0.2* (t.steps[ijump]-t.steps[0])

    visualize.scene()
    visualize.box(t.read(0).cell)
    for i0 in t.samples[0::ijump]:
        if i0 + i >= len(t.samples):
            break
        s0f = ti0.read(i0)
        s0 = t0.read(i0)
        s1 = t.read(i0+i)

        dr = []
        for p0f, p0, p1 in zip(s0f.particle, s0.particle, s1.particle):            
            dr.append(sqrt(sum((p1.position-p0.position)**2)))
        visualize.show(s0f.particle, dr, rscale=0.35, rcut=0.0)
        if not out_dir is None:
            visualize.dump(k, out_dir)
        k+=1
    if movie:
        visualize.movie(out_dir)

if __name__ == '__main__':
    #angle(sys.argv[1], float(sys.argv[2]))
    displacement(trj.Trajectory(sys.argv[1]), int(sys.argv[2]), '/tmp/visualize')
