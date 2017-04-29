#!/usr/bin/env python

import sys
import os
import unittest
import numpy
from atooms import trajectory
import postprocessing
from postprocessing.helpers import filter_species

def filter_id(t, s, id):
    # TODO: should we do a copy or modify the system in place? Looks like modify it's ok
    nop = [p for p in s.particle if p.id != id]
    for n in nop:
        s.particle.remove(n)
    return s

def filter_random(t, s, n):
    """Keep only n particles"""
    import random
    ids = random.sample(xrange(len(s.particle)), len(s.particle)-n)
    nop = [s.particle[i] for i in ids]
    for p in nop:
        s.particle.remove(p)
    return s

def filter_selected_ids(t, s, ids):
    """Keep only selected ids of particles"""
    nop = [s.particle[i] for i in ids]
    for p in nop:
        s.particle.remove(p)
    return s

def deviation(x, y):
    return (numpy.sum((x-y)**2)/len(x))**0.5

class TestRealSpace(unittest.TestCase):

    def setUp(self):
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')

    def test_msd_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.Sliced(trajectory.TrajectoryHDF5(f), slice(0,1000,1))
        for i in [1, 2]:
            t = trajectory.Filter(ts, filter_id, i)
            p = postprocessing.MeanSquareDisplacement(t, [0.0, 3.0, 45.0, 90])
            p.compute()            
            ref_grid = numpy.array([0, 3.0, 45.0, 90.0])
            if i == 1:
                ref_value = numpy.array([0.0, 0.126669, 1.21207, 2.16563])
            else:
                ref_value = numpy.array([0.0, 0.220299, 2.31111, 4.37561])
            self.assertLess(deviation(p.grid, ref_grid), 4e-2)
            self.assertLess(deviation(p.value, ref_value), 4e-2)

    def test_msd_partial_filter(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.Sliced(trajectory.TrajectoryHDF5(f), slice(0,1000,1))
        for i in [1, 2]:
            p = postprocessing.MeanSquareDisplacement(ts, [0.0, 3.0, 45.0, 90])
            p.add_filter(filter_species, i)
            p.compute()
            ref_grid = numpy.array([0, 3.0, 45.0, 90.0])
            if i == 1:
                ref_value = numpy.array([0.0, 0.126669, 1.21207, 2.16563])
            else:
                ref_value = numpy.array([0.0, 0.220299, 2.31111, 4.37561])
            self.assertLess(deviation(p.grid, ref_grid), 4e-2)
            self.assertLess(deviation(p.value, ref_value), 4e-2)

    def test_gr_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        ref = {}
        ref[(1,1)] = numpy.array([ 0.,          0.00675382,  0.27087136,  1.51486318])
        ref[(2,2)] = numpy.array([ 0.31065645,  0.51329066,  0.67485665,  0.78039485])
        ref[(1,2)] = numpy.array([ 4.25950671,  3.86572027,  2.70020052,  1.78935426])
        for i in [1, 2]:
            p = postprocessing.RadialDistributionFunction(ts)
            p.add_filter(filter_species, i)
            r, gr = p.compute()
            self.assertLess(deviation(gr[21:25], ref[(i,i)]), 4e-2)

        p = postprocessing.RadialDistributionFunction(ts)
        p.add_filter(filter_species, 1)
        p.add_filter(filter_species, 2)
        r, gr = p.compute()
        self.assertLess(deviation(gr[21:25], ref[(1,2)]), 4e-2)

    def test_gr_partial_2(self):
        from postprocessing.partial import Partial
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        ref = {}
        ref[(1,1)] = numpy.array([ 0.,          0.00675382,  0.27087136,  1.51486318])
        ref[(2,2)] = numpy.array([ 0.31065645,  0.51329066,  0.67485665,  0.78039485])
        ref[(1,2)] = numpy.array([ 4.25950671,  3.86572027,  2.70020052,  1.78935426])
        
        gr = Partial(postprocessing.RadialDistributionFunction, [1,2], ts)
        gr.compute()
        for ab in [(1,1), (1,2), (2,2)]:
            self.assertLess(deviation(gr.partial[ab].value[21:25], ref[ab]), 4e-2)

class TestFourierSpace(unittest.TestCase):

    def setUp(self):
        import random
        random.seed(10)
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')

    def test_sk(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        t = trajectory.TrajectoryHDF5(f)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.compute()
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.04)

    def test_sk_fixgrid(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        t = trajectory.TrajectoryHDF5(f)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10])
        p.compute()
        ref_value = numpy.array([0.083411717745282138, 2.76534619194135, 0.67129958432631986])
        self.assertLess(deviation(p.value, ref_value), 0.04)

    def test_sk_variable_cell(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        p = postprocessing.StructureFactor(ts, range(1,10))
        p.compute()
        t = trajectory.AffineDeformation(ts, 1e-3)
        p1 = postprocessing.StructureFactor(t, range(1,10))
        p1.compute()
        # for x, y, z, w in zip(p.grid, p1.grid, p.value, p1.value):
        #     print x, y, z, w

    def test_sk_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        for i in [1, 2]:
            t = trajectory.Filter(ts, filter_id, i)
            p = postprocessing.StructureFactor(t, [4, 7.3, 10])
            p.compute()
            if i == 1:
                ref_value = numpy.array([0.078218154990600072, 2.8964366727431656, 0.54336352118808429])
            else:
                ref_value = numpy.array([0.86716496871363735, 0.86986885176760842, 0.98112175463699136])
            self.assertLess(deviation(p.value, ref_value), 1e-2)

    def test_sk_random(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        t = trajectory.Filter(ts, filter_random, 75)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10, 30.0], nk=40)
        p.compute()

    @unittest.skip('Broken test')
    def test_fkt_random(self):
        import random
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        ts = trajectory.TrajectoryHDF5(f)
        s = ts[0]
        ids = random.sample(xrange(len(s.particle)), len(s.particle))
        t = trajectory.Filter(ts, filter_selected_ids, ids)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40)
        p.compute()
        ps = postprocessing.IntermediateScattering(ts, [4, 7.3, 10], nk=40)
        ps.compute()

    def test_fkt_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.h5')
        t = trajectory.TrajectoryHDF5(f)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40)
        p.add_filter(filter_species, 1)
        p.compute()

if __name__ == '__main__':
    unittest.main()


