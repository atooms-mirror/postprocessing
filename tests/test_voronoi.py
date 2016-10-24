#!/usr/bin/env python

import unittest
import math
from atooms.trajectory import Trajectory
from atooms.plugins import voronoi

class TestVoronoi(unittest.TestCase):

    def setUp(self):
        pass

    def test_read(self):
        t = voronoi.TrajectoryVoronoi('reference/wahn_voronoi.xyz')
        v0 = []
        for i in t.samples:
            s = t.read_sample(i)
            v = s.voronoi
            v0.append(v[0].signature)
        ref = [(0, 1, 10, 4), (0, 3, 6, 7), (1, 1, 8, 3, 1), (0, 2, 8, 4), (1, 5, 2, 5, 3), (0, 1, 10, 4), (0, 2, 8, 4), (0, 3, 7, 5, 1), (2, 1, 7, 4, 1, 1), (0, 1, 10, 4), (0, 3, 6, 4), (0, 2, 8, 5), (0, 2, 8, 6), (2, 3, 3, 5, 1, 1), (0, 2, 8, 5), (0, 2, 8, 5), (0, 1, 10, 4), (4, 3, 2, 3, 2, 3), (0, 1, 10, 3), (1, 1, 8, 3, 1)]
        self.assertEqual(ref, v0)

    def test_voro(self):
        """Test of polydisperse hard spheres"""
        finp = '/home/coslo/projects/polydisperse_swap/data/misaki/const_volume/EQ/N1000/phi0.643/conv-configs/config.h5'
        fout = '/tmp/voro.xyz'
        t = Trajectory(finp)
        s = t[0]
        vt = voronoi.VoronoiTessellation(s, fout)
        # Sum of Voronoi volumes should equal V
        dV = sum([v.volume for v in vt.polyhedra]) - s.cell.volume
        self.assertLess(abs(dV) / s.cell.volume, 1e-6)
        # Accumulated free volume should give back 1-phi (packing fraction)
        phi = sum([4./3*math.pi*p.radius**3 for p in s.particle]) / s.cell.volume
        free_volume = sum([v.volume - 4/3.*math.pi*p.radius**3 for p, v in zip(s.particle, vt.polyhedra)])
        one = free_volume / s.cell.volume + phi
        self.assertLess(abs(one-1), 1e-6)
            
if __name__ == '__main__':
    unittest.main(verbosity=0)


