#!/usr/bin/env python

import sys
import os
import random
import unittest
import numpy
from atooms import trajectory
import atooms.postprocessing as postprocessing
from atooms.postprocessing.helpers import filter_species

# def filter_id(s, id):
#     # TODO: should we do a copy or modify the system in place? Looks like modify it's ok
#     nop = [p for p in s.particle if p.species != id]
#     for n in nop:
#         s.particle.remove(n)
#     return s

def filter_random(s, n):
    """Keep only n particles"""
    ids = random.sample(range(len(s.particle)), len(s.particle)-n)
    nop = [s.particle[i] for i in ids]
    for p in nop:
        s.particle.remove(p)
    return s

def filter_selected_ids(s, ids):
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
        ref_grid = numpy.array([0, 3.0, 45.0, 90.0])
        ref_value = {'A': numpy.array([0.0, 0.126669, 1.21207, 2.16563]),
                     'B': numpy.array([0.0, 0.220299, 2.31111, 4.37561])}
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        for i in ['A', 'B']:
            with trajectory.Sliced(trajectory.TrajectoryXYZ(f), slice(0, 1000, 1)) as t:
                t.add_callback(filter_species, i)
                p = postprocessing.MeanSquareDisplacement(t, [0.0, 3.0, 45.0, 90])
                import warnings
                warnings.simplefilter('ignore', RuntimeWarning)
                p.compute()
                self.assertLess(deviation(p.grid, ref_grid), 4e-2)
                self.assertLess(deviation(p.value, ref_value[i]), 4e-2)

    def test_msd_partial_filter(self):
        ref_grid = numpy.array([0, 3.0, 45.0, 90.0])
        ref_value = {'A': numpy.array([0.0, 0.126669, 1.21207, 2.16563]),
                     'B': numpy.array([0.0, 0.220299, 2.31111, 4.37561])}
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.Sliced(trajectory.TrajectoryXYZ(f), slice(0, 1000, 1))
        for i in ['A', 'B']:
            p = postprocessing.MeanSquareDisplacement(ts, [0.0, 3.0, 45.0, 90])
            p.add_filter(filter_species, i)
            p.compute()
            self.assertLess(deviation(p.grid, ref_grid), 4e-2)
            self.assertLess(deviation(p.value, ref_value[i]), 4e-2)
        ts.close()

    def test_gr_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        ref = {}
        ref[('A', 'A')] = numpy.array([ 0.,          0.00675382,  0.27087136,  1.51486318])
        ref[('B', 'B')] = numpy.array([ 0.31065645,  0.51329066,  0.67485665,  0.78039485])
        ref[('A', 'B')] = numpy.array([ 4.25950671,  3.86572027,  2.70020052,  1.78935426])
        for i in ['A', 'B']:
            p = postprocessing.RadialDistributionFunction(ts)
            p.add_filter(filter_species, i)
            r, gr = p.compute()
            self.assertLess(deviation(gr[21:25], ref[(i, i)]), 4e-2)

        p = postprocessing.RadialDistributionFunction(ts)
        p.add_filter(filter_species, 'A')
        p.add_filter(filter_species, 'B')
        r, gr = p.compute()
        self.assertLess(deviation(gr[21:25], ref[('A', 'B')]), 4e-2)

    def test_gr_partial_2(self):
        from atooms.postprocessing.partial import Partial
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        ref = {}
        ref[('A', 'A')] = numpy.array([ 0.,          0.00675382,  0.27087136,  1.51486318])
        ref[('B', 'B')] = numpy.array([ 0.31065645,  0.51329066,  0.67485665,  0.78039485])
        ref[('A', 'B')] = numpy.array([ 4.25950671,  3.86572027,  2.70020052,  1.78935426])
        
        gr = Partial(postprocessing.RadialDistributionFunction, ['A', 'B'], ts)
        gr.compute()
        for ab in [('A', 'A'), ('A', 'B'), ('B', 'B')]:
            self.assertLess(deviation(gr.partial[ab].value[21:25], ref[ab]), 4e-2)

class TestFourierSpace(unittest.TestCase):

    def setUp(self):
        random.seed(10)
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')

    def test_sk(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.compute()
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.04)

    def test_sk_opti(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactorOptimized(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.compute()
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.05)

    def test_sk_update(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.do(update=False)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.do(update=True)
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.04)

    def test_sk_fixgrid(self):
        # TODO: this test fails with python 3 (small deviations)
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10])
        p.compute()
        ref_value = numpy.array([0.083411717745282138, 2.76534619194135, 0.67129958432631986])
        self.assertLess(deviation(p.value, ref_value), 0.04)

    def test_sk_variable_cell(self):
        # TODO: this test has no assertion
        def deformation(s, scale=0.01):
            # Note this random scaling changes every time read is called,
            # even for the same sample
            x = 1 + (random.random()-0.5) * scale
            s.cell.side *= x
            for p in s.particle:
                p.position *= x
            return s
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        with trajectory.TrajectoryXYZ(f) as t:
            p = postprocessing.StructureFactor(t, list(range(1,10)))
            p.compute()
        with trajectory.TrajectoryXYZ(f) as t:
            t.add_callback(deformation, 1e-3)
            p = postprocessing.StructureFactor(t, list(range(1,10)))
            p.compute()

    def test_sk_partial(self):
        # TODO: this test fails with python 3 (small deviations)
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ref_value = {'A': numpy.array([0.078218, 2.896436, 0.543363]),
                     'B': numpy.array([0.867164, 0.869868, 0.981121]),
                     'AB': numpy.array([-0.1907, 0.399360, 0.050480])}
        for species in ['A', 'B']:
            with trajectory.TrajectoryXYZ(f) as t:
                t.add_callback(filter_species, species)
                p = postprocessing.StructureFactor(t, [4, 7.3, 10])
                p.compute()
                self.assertLess(deviation(p.value, ref_value[species]), 1e-2)

        with trajectory.TrajectoryXYZ(f) as t:
            sk = postprocessing.Partial(postprocessing.StructureFactor, ['A', 'B'], t, [4, 7.3, 10])
            sk.compute()
            self.assertLess(deviation(sk.partial[('A', 'A')].value, ref_value['A']), 1e-2)
            self.assertLess(deviation(sk.partial[('B', 'B')].value, ref_value['B']), 1e-2)
            self.assertLess(deviation(sk.partial[('A', 'B')].value, ref_value['AB']), 1e-2)

    def test_sk_random(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_random, 75)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10, 30.0], nk=40)
        p.compute()

    def test_sk_field(self):
        """
        Test that S(k) with a field that is 0 if id=A and 1 if id=B gives
        the BB partial structure factor.
        """
        # TODO: this test fails with python 3 because of a weird issue with xyz trajectory in atooms (_fallback)
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ff = os.path.join(self.reference_path, 'kalj-small-field.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10], trajectory_field=ff)
        p.compute()
        # We multiply by x because the S(k) is normalized to 1/N
        from atooms.system.particle import composition
        x = composition(t[0].particle)['B'] / float(len(t[0].particle))
        ref_value = x * numpy.array([0.86716496871363735, 0.86986885176760842, 0.98112175463699136])
        self.assertLess(deviation(p.value, ref_value), 1e-2)

    @unittest.skip('Broken test')
    def test_fkt_random(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        with trajectory.TrajectoryXYZ(f) as t:
            s = t[0]
            ids = random.sample(range(len(s.particle)), len(s.particle))
            t.add_callback(filter_selected_ids, ids)
            p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40)
            p.compute()

    def test_fkt_partial(self):
        # TODO: add check
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40)
        p.add_filter(filter_species, 'A')
        p.compute()

    def test_fskt_partial(self):
        # TODO: add check
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.SelfIntermediateScattering(t, [4, 7.3, 10], nk=40, norigins=0.2)
        p.add_filter(filter_species, 'A')
        p.compute()

if __name__ == '__main__':
    unittest.main()


