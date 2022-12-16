#!/usr/bin/env python

import sys
import os
import random
import unittest
import numpy
from atooms import trajectory
import atooms.postprocessing as postprocessing
from atooms.postprocessing.helpers import filter_species


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

def filter_2d(s):
    s.cell.side = s.cell.side[0:2]
    s.cell.center = s.cell.center[0:2]
    for p in s.particle:
        p.position = p.position[0:2]
        p.velocity = p.velocity[0:2]
    return s

class Test(unittest.TestCase):

    def test_name(self):
        reference_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
        default = postprocessing.core.pp_output_path
        postprocessing.core.pp_output_path = '{trajectory.filename}.pp.{short_name}.{tag_description}'
        corr = postprocessing.SelfIntermediateScattering(os.path.join(reference_path, 'trajectory.xyz'))
        self.assertEqual(os.path.basename(corr._output_file), 'trajectory.xyz.pp.F_s(k,t).the_whole_system')
        self.assertEqual(corr.grid_name, ['k', 't'])
        postprocessing.core.pp_output_path = default
        corr.trajectory.close()

class TestRealSpace(unittest.TestCase):

    def setUp(self):
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')

    def test_msd_partial(self):
        import warnings
        warnings.simplefilter('ignore', RuntimeWarning)
        from atooms.postprocessing.partial import Partial
        ref_grid = numpy.array([0, 3.0, 45.0, 90.0])
        ref_value = {'A': numpy.array([0.0, 0.12708190185520277, 1.2131017102523827, 2.1696872045992075]),
                     'B': numpy.array([0.0, 0.22107274740436153, 2.3018226393609473, 4.354207960272026])}
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        for i in ['A', 'B']:
            with trajectory.Sliced(trajectory.TrajectoryXYZ(f), slice(0, 1000, 1)) as t:
                p = Partial(postprocessing.MeanSquareDisplacement, ['A', 'B'], t, [0.0, 3.0, 45.0, 90], fix_cm=True)
                p.compute()
                self.assertLess(deviation(p.partial[i].grid, ref_grid), 4e-2)
                self.assertLess(deviation(p.partial[i].value, ref_value[i]), 4e-2)

        f = os.path.join(self.reference_path, 'kalj-small-unfolded.xyz')
        for i in ['A', 'B']:
            with trajectory.Sliced(trajectory.TrajectoryXYZ(f), slice(0, 1000, 1)) as t:
                p = postprocessing.MeanSquareDisplacement(t, [0.0, 3.0, 45.0, 90], fix_cm=True)
                p.add_filter(filter_species, i)
                p.compute()
                self.assertLess(deviation(p.grid, ref_grid), 4e-2)
                self.assertLess(deviation(p.value, ref_value[i]), 4e-2)

        # This calculation is correct, because filter is applied after unfolding
        # If we apply the filter to the trajectory, there are slight differences for small particles
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.Sliced(trajectory.TrajectoryXYZ(f), slice(0, 1000, 1))
        for i in ['A', 'B']:
            p = postprocessing.MeanSquareDisplacement(ts, [0.0, 3.0, 45.0, 90], fix_cm=True)
            p.add_filter(filter_species, i)
            p.compute()
            self.assertLess(deviation(p.grid, ref_grid), 4e-2)
            self.assertLess(deviation(p.value, ref_value[i]), 4e-2)
        ts.close()

    def _test_gr_partial(self, cls):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        ref = {}
        ref[('A', 'A')] = numpy.array([0., 0.00675382, 0.27087136, 1.51486318])
        ref[('B', 'B')] = numpy.array([0.31065645, 0.51329066, 0.67485665, 0.78039485])
        ref[('A', 'B')] = numpy.array([4.25950671, 3.86572027, 2.70020052, 1.78935426])
        for i in ['A', 'B']:
            p = cls(ts)
            p.add_filter(filter_species, i)
            r, gr = p.compute()
            self.assertLess(deviation(gr[21:25], ref[(i, i)]), 4e-2)

        p = cls(ts)
        p.add_filter(filter_species, 'A')
        p.add_filter(filter_species, 'B')
        r, gr = p.compute()
        self.assertLess(deviation(gr[21:25], ref[('A', 'B')]), 4e-2)
        ts.close()

    def test_gr_partial_big(self):
        f = os.path.join(self.reference_path, 'ka_N20000.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        isp = 'B'
        res = {}
        grid = numpy.linspace(0.1, 5.0, int(5.0/0.1))
        p = postprocessing.RadialDistributionFunctionFast(ts, dr=0.1, rgrid=grid, rmax=-1)
        p.add_filter(filter_species, isp)
        res[-1] = p.compute()
        # TODO: fix issue with arbitrary grid
        p1 = postprocessing.RadialDistributionFunctionFast(ts, rmax=5)
        p1.add_filter(filter_species, isp)
        res[5] = p1.compute()
        # p.show(now=False)
        # p1.show()
        #self.assertLess(deviation(res[-1][1], res[5][1]), 1e-3)
        ts.close()

    def test_gr_big(self):
        f = os.path.join(self.reference_path, 'ka_N20000.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        p = postprocessing.RadialDistributionFunctionFast(ts, rmax=5)
        p.compute()
        #self.assertLess(deviation(res[-1][1], res[5][1]), 1e-6)
        ts.close()

    def test_gr_partial_fast(self):
        # This will test fast if available
        self._test_gr_partial(postprocessing.RadialDistributionFunction)

    def test_gr_partial_legacy(self):
        self._test_gr_partial(postprocessing.RadialDistributionFunctionLegacy)

    def test_gr_partial_2(self):
        from atooms.postprocessing.partial import Partial
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ref = {}
        ref[('A', 'A')] = numpy.array([0., 0.00675382, 0.27087136, 1.51486318])
        ref[('B', 'B')] = numpy.array([0.31065645, 0.51329066, 0.67485665, 0.78039485])
        ref[('A', 'B')] = numpy.array([4.25950671, 3.86572027, 2.70020052, 1.78935426])
        with trajectory.TrajectoryXYZ(f) as ts:
            gr = Partial(postprocessing.RadialDistributionFunction, ['A', 'B'], ts)
            gr.compute()
            for ab in [('A', 'A'), ('A', 'B'), ('B', 'B')]:
                self.assertLess(deviation(gr.partial[ab].value[21:25], ref[ab]), 4e-2)

    def test_gr_filter(self):
        from atooms.postprocessing.filter import Filter
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        ref = {}
        ref[('A', 'A')] = numpy.array([0., 0.00675382, 0.27087136, 1.51486318])
        ref[('B', 'B')] = numpy.array([0.31065645, 0.51329066, 0.67485665, 0.78039485])
        ref[('A', 'B')] = numpy.array([4.25950671, 3.86572027, 2.70020052, 1.78935426])

        gr = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "A", species == "A"')
        gr.compute()
        self.assertLess(deviation(gr.value[21:25], ref[('A', 'A')]), 4e-2)

        gr = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "A", species == "B"')
        gr.compute()
        self.assertLess(deviation(gr.value[21:25], ref[('A', 'B')]), 4e-2)

        gr = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "B", species == "B"')
        gr.compute()
        self.assertLess(deviation(gr.value[21:25], ref[('B', 'B')]), 4e-2)
        ts.close()

    def test_gr_filter_2(self):
        from atooms.postprocessing.filter import Filter
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ts = trajectory.TrajectoryXYZ(f)

        gr_AX = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "A"')
        gr_AX.compute()

        gr_AA = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "A", species == "A"')
        gr_AA.compute()

        gr_AB = Filter(postprocessing.RadialDistributionFunction(ts), 'species == "A", species == "B"')
        gr_AB.compute()

        self.assertLess(deviation(gr_AX.value[15:25], 0.8 * gr_AA.value[15:25] + 0.2 * gr_AB.value[15:25]), 0.04)
        ts.close()

    def test_gr_partial_2d_open(self):
        from atooms.postprocessing.partial import Partial
        import numpy
        def fix_2d(system):
            cm = system.cm_position[: 2]
            system.cell.periodic = numpy.array([False, False])
            system.cell.side = system.cell.side[0: 2]
            for p in system.particle:
                p.position = p.position[0: 2] - cm
            return system
        f = os.path.join('atoms_coordinates_test.xyz')
        ts = trajectory.TrajectoryXYZ(f)
        ts.add_callback(fix_2d)        
        gr = Partial(postprocessing.RadialDistributionFunction, ['0', '1'], ts, rmax=2.0)
        gr.compute()

        import matplotlib.pyplot as plt
        plt.plot(gr.partial[('0', '0')].grid[1:], gr.partial[('0', '0')].value[1:])
        plt.plot(gr.partial[('0', '1')].grid[1:], gr.partial[('0', '1')].value[1:])
        plt.plot(gr.partial[('1', '1')].grid[1:], gr.partial[('1', '1')].value[1:])
        plt.grid()
        # p = ts[0].dump('pos')
        # plt.plot(p[:, 0], p[:, 1], 'o')
        plt.show()
        ts.close()

        
class TestFourierSpace(unittest.TestCase):

    def setUp(self):
        random.seed(10)
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')

    def _test_sk(self, cls):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = cls(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.compute()
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.04)
        # Trivial report dump (no check)
        txt = p.report()
        t.close()

    def test_sk(self):
        self._test_sk(postprocessing.StructureFactor)

    def test_sk_fast(self):
        try:
            import atooms.postprocessing.fourierspace_wrap
            self._test_sk(postprocessing.StructureFactorFast)
        except ImportError:
            self.skipTest('missing f90')

    def test_sk_update(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.do(update=False)
        p = postprocessing.StructureFactor(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.do(update=True)
        ref_value = numpy.array([0.075820086512828039, 0.065300213310725302, 0.082485082309989494])
        self.assertLess(deviation(p.value, ref_value), 0.04)
        t.close()

    def test_sk_fixgrid(self):
        # TODO: this test fails with python 3 (small deviations)
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10])
        p.compute()
        ref_value = numpy.array([0.083411717745282138, 2.76534619194135, 0.67129958432631986])
        self.assertLess(deviation(p.value, ref_value), 0.08)
        t.close()

    def test_kvectors(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10], nk=10)
        p.compute()
        q = postprocessing.StructureFactor(t)
        q.kvectors = p.kvectors
        self.assertEqual(q.kvectors, p.kvectors)
        self.assertEqual(q.kgrid, p.kgrid)
        q.compute()
        self.assertEqual(q.grid, p.grid)
        self.assertEqual(q.value, p.value)
        t.close()

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
            p = postprocessing.StructureFactor(t, [4, 7.3, 10])
            p.compute()
        with trajectory.TrajectoryXYZ(f) as t:
            t.add_callback(deformation, 1e-3)
            p = postprocessing.StructureFactor(t, [4, 7.3, 10])
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
                self.assertLess(deviation(p.value, ref_value[species]), 1e-1)

        with trajectory.TrajectoryXYZ(f) as t:
            sk = postprocessing.Partial(postprocessing.StructureFactor, ['A', 'B'], t, [4, 7.3, 10])
            sk.compute()
            self.assertLess(deviation(sk.partial[('A', 'A')].value, ref_value['A']), 1e-1)
            self.assertLess(deviation(sk.partial[('B', 'B')].value, ref_value['B']), 1e-1)
            self.assertLess(deviation(sk.partial[('A', 'B')].value, ref_value['AB']), 1e-1)

    def test_sk_random(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_random, 75)
        p = postprocessing.StructureFactor(t, [4, 7.3, 10, 30.0], nk=40)
        p.compute()
        t.close()

    def test_sk_field(self):
        """
        Test that S(k) with a field that is 0 if id=A and 1 if id=B gives
        the BB partial structure factor.
        """
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ff = os.path.join(self.reference_path, 'kalj-small-field.xyz')
        th = trajectory.TrajectoryXYZ(f)
        tt = trajectory.TrajectoryXYZ(ff)
        p = postprocessing.StructureFactor(th, [4, 7.3, 10])
        p.add_weight(trajectory=tt, field='field_B')
        p.compute()
        # We multiply by x because the S(k) is normalized to 1/N
        from atooms.system.particle import composition
        x = composition(th[0].particle)['B'] / float(len(th[0].particle))
        ref_value = x * numpy.array([0.86716496871363735, 0.86986885176760842, 0.98112175463699136])
        self.assertLess(deviation(p.value, ref_value), 1e-2)
        th.close()
        tt.close()

    def test_sk_field_partial(self):
        """
        Test that weight works with partial correlation
        """
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        ff = os.path.join(self.reference_path, 'kalj-small-field.xyz')
        th = trajectory.TrajectoryXYZ(f)
        p = postprocessing.Partial(postprocessing.StructureFactor, ['A', 'B'], th, [4, 7.3, 10])
        from atooms.postprocessing.helpers import copy_field
        from atooms.trajectory import TrajectoryXYZ
        p.add_weight(trajectory=trajectory.TrajectoryXYZ(ff), field='field_B')
        p.compute()
        from atooms.system.particle import composition
        ref_value = numpy.array([0.86716496871363735, 0.86986885176760842, 0.98112175463699136])
        zeros = numpy.zeros(3)
        self.assertLess(deviation(p.partial[('B', 'B')].value, ref_value), 2.3e-2)
        self.assertLess(deviation(p.partial[('A', 'A')].value, zeros), 2e-2)
        th.close()

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
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40)
        p.add_filter(filter_species, 'A')
        p.compute()
        p.analyze()
        tau = []
        for key in sorted(p.analysis['relaxation times tau']):
            tau.append(p.analysis['relaxation times tau'][key])
        self.assertLess(abs(tau[0] - 2.2792074711157104), 0.4)
        self.assertLess(abs(tau[1] - 5.8463508731564975), 0.4)
        self.assertLess(abs(tau[2] - 0.85719855804743605), 0.4)
        t.close()

    def test_fkt_nonorm_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=40, normalize=False)
        p.add_filter(filter_species, 'A')
        p.compute()
        p.analyze()
        tau = []
        for key in sorted(p.analysis['relaxation times tau']):
            tau.append(p.analysis['relaxation times tau'][key])
        self.assertLess(abs(tau[0] - 2.2792074711157104), 0.4)
        self.assertLess(abs(tau[1] - 5.8463508731564975), 0.4)
        self.assertLess(abs(tau[2] - 0.85719855804743605), 0.4)
        t.close()

    def test_fskt_partial(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        for cls in [postprocessing.SelfIntermediateScatteringLegacy, postprocessing.SelfIntermediateScatteringFast]:
            p = cls(t, [4, 7.3, 10], nk=40, norigins=0.2)
            p.add_filter(filter_species, 'A')
            p.compute()
            p.analyze()
            tau = []
            for key in sorted(p.analysis['relaxation times tau']):
                tau.append(p.analysis['relaxation times tau'][key])
            self.assertLess(abs(tau[0] - 14.081572329287619), 0.04)
            self.assertLess(abs(tau[1] - 3.1034088042905967), 0.04)
            self.assertLess(abs(tau[2] - 0.97005294966138289), 0.04)
        t.close()

    def test_chi4_overlap(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        with trajectory.TrajectoryXYZ(f) as th:
            tgrid = postprocessing.helpers.logx_grid(0.0, th.total_time * 0.5, 10)
            fct = postprocessing.Susceptibility(postprocessing.SelfOverlap, th, tgrid=tgrid)
            fct.compute()
            ref = postprocessing.Chi4SelfOverlap(th, tgrid=tgrid)
            ref.compute()
            self.assertLess(deviation(numpy.array(ref.value), numpy.array(fct.value)), 0.1)

    def test_fskt_2d(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_2d)
        p = postprocessing.SelfIntermediateScattering(t, [4, 7.3, 10], nk=10)
        p.add_filter(filter_species, 'A')
        p.compute()
        p.analyze()
        tau = []
        for key in sorted(p.analysis['relaxation times tau']):
            tau.append(p.analysis['relaxation times tau'][key])
        self.assertLess(abs(tau[0] - 13.48342847723456), 0.04)
        self.assertLess(abs(tau[1] - 3.07899513664358), 0.04)
        self.assertLess(abs(tau[2] - 0.9802163934982774), 0.04)
        t.close()

    def test_fskt_legacy_2d(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_2d)
        p = postprocessing.SelfIntermediateScatteringLegacy(t, [4, 7.3, 10], nk=10)
        p.add_filter(filter_species, 'A')
        p.compute()
        p.analyze()
        tau = []
        for key in sorted(p.analysis['relaxation times tau']):
            tau.append(p.analysis['relaxation times tau'][key])
        self.assertLess(abs(tau[0] - 13.48342847723456), 0.04)
        self.assertLess(abs(tau[1] - 3.07899513664358), 0.04)
        self.assertLess(abs(tau[2] - 0.9802163934982774), 0.04)
        t.close()

    def test_fkt_2d(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_2d)
        p = postprocessing.IntermediateScattering(t, [4, 7.3, 10], nk=100)
        p.add_filter(filter_species, 'A')
        p.compute()
        p.analyze()
        tau = []
        for key in sorted(p.analysis['relaxation times tau']):
            tau.append(p.analysis['relaxation times tau'][key])
        self.assertLess(abs(tau[0] - 1.1341521365187757), 0.04)
        self.assertLess(abs(tau[1] - 5.83114954720099), 0.04)
        self.assertLess(abs(tau[2] - 0.859950963462569), 0.04)
        t.close()

    def test_sk_2d(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        t.add_callback(filter_2d)
        p = postprocessing.StructureFactorLegacy(t, kmin=-1, kmax=4, ksamples=3, dk=0.2)
        p.compute()
        ref_value = numpy.array([[0.06899986228704291, 0.0629709003150001, 0.07397620251792263]])
        self.assertLess(deviation(p.value, ref_value), 0.04)
        t.close()

    def test_fourierspace_kgrid_unsorted(self):
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p_sorted = postprocessing.SelfIntermediateScatteringLegacy(t, [4, 7.3, 10], nk=40, norigins=0.2)
        p_unsorted = postprocessing.SelfIntermediateScatteringLegacy(t, [10, 4, 7.3], nk=40, norigins=0.2)

        p_sorted.compute()
        p_sorted.analyze()
        p_unsorted.compute()
        p_unsorted.analyze()

        self.assertEqual(p_sorted.kgrid, p_unsorted.kgrid)
        self.assertLess(deviation(numpy.array(p_sorted.value), numpy.array(p_unsorted.value)), 1e-14)
        t.close()

    def test_ba(self):
        """Check on cut-off for neighbors determination"""
        f = os.path.join(self.reference_path, 'kalj-small.xyz')
        t = trajectory.TrajectoryXYZ(f)
        p = postprocessing.BondAngleDistribution(t)
        p.compute()
        q = postprocessing.BondAngleDistribution(t, rcut=p.rcut)
        q.compute()
        self.assertLess(deviation(p.value, q.value), 0.001)
        t.close()

    @unittest.skip('Broken test')
    def test_gr_crop(self):
        # TODO: fix this test
        import numpy
        import atooms.postprocessing as pp
        import atooms.trajectory as trj
        import atooms.system
        numpy.random.seed(1)
        N = 5000
        L = 1000.0
        pos = numpy.random.random([N, 3]) * L
        system = atooms.system.System()
        system.particle = [atooms.system.Particle(position=pos[i, :]) for i in range(N)]

        def center(system):
            cm = system.cm_position
            for p in system.particle:
                p.position[:] -= cm
            return system

        def bounding_box(system):
            L = []
            periodic = numpy.ndarray(system.number_of_dimensions, dtype=bool)
            periodic[:] = False
            for axis in range(system.number_of_dimensions):
                L.append(1.01 * 2 * numpy.max([abs(p.position[axis]) for p in system.particle]))
            system.cell = atooms.system.Cell(L, periodic=periodic)
            return system

        def crop(system, L):
            new = []
            for p in system.particle:
                if numpy.all(numpy.abs(p.position) < L / 2):
                    new.append(p)
            system.particle = new
            return system

        rmax = L / 10
        th = trj.TrajectoryRam()
        th[0] = system
        th.add_callback(center)
        th.add_callback(crop, numpy.array([L, L, L]))
        th.add_callback(bounding_box)

        gr = pp.RadialDistributionFunction(th, dr=1.0, rmax=rmax)
        gr.compute()
        self.assertLess(abs(gr.value[-1] - 1.0), 0.1)


if __name__ == '__main__':
    unittest.main()
