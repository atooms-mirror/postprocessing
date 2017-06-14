#!/usr/bin/env python

import unittest
import numpy
from atooms.trajectory import TrajectoryNeighbors, Trajectory
from postprocessing.boo import BondOrientationalOrder, periodic_vector


class TestBOO(unittest.TestCase):

    def setUp(self):
        pass

    def _load(self, structure):
        if structure == 'bcc':
            self.t = Trajectory('data/lj_bcc.xyz')
            self.tn = TrajectoryNeighbors('data/lj_bcc.xyz.voronoi.xyz.neigh')
            self.ref = {6: 0.510688230857, 4: 0.0363696483727, 8: 0.429322472922}
        elif structure == 'fcc':
            self.t = Trajectory('data/lj_fcc.xyz.min')
            self.tn = TrajectoryNeighbors('data/lj_fcc.xyz.min.voronoi.xyz.neigh')
            self.ref = {6: 0.574524259714, 4: 0.190940653956, 8: 0.403914561085}
        elif structure == 'fluid':
            self.t = Trajectory('data/lj.xyz')
            self.tn = TrajectoryNeighbors('data/lj.xyz.neigh')
            #self.ref = {6: 0.2731} # voronoi
            self.ref = {6: 0.33194535071671305}

    def _test(self):
        s = self.t[0]
        np = len(s.particle)
        box = s.cell.side
        n = self.tn[0]
        boo = BondOrientationalOrder(s.particle, n.neighbors, box)
        for l in self.ref:
            self.assertAlmostEqual(self.ref[l], numpy.average(boo.ql(l)))

    # def test_fcc(self):
    #     self._load('fcc')
    #     self._test()

    def test_bcc(self):
        self._load('bcc')
        self._test()

    def test_fluid(self):
        self._load('fluid')
        self._test()
        

if __name__ == '__main__':
    unittest.main(verbosity=0)


