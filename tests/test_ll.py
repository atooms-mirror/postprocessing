#!/usr/bin/env python

import sys
import os
import random
import unittest
import numpy
from atooms.postprocessing import linkedcells


def deviation(x, y):
    return (numpy.sum((x-y)**2)/len(x))**0.5

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_basic(self):
        N = 4
        ndim = 3
        L = 6.0
        rcut = 1.0
        side = numpy.ndarray(ndim)
        pos = numpy.ndarray((N, ndim))
        side[:] = [L, L, L]
        pos[0, :] = [-0.5, 0.0, 0.0]
        pos[1, :] = [0.5, 0.0, 0.0]
        pos[2, :] = [1.5, 0.0, 0.0]
        pos[3, :] = [1.5, 1.5, 0.0]
        lc = linkedcells.LinkedCells(rcut=rcut)
        nn, num = lc.compute(side, pos, as_array=True, newton=False)
        self.assertEqual(set(nn[0][0: num[0]]), set([1]))
        self.assertEqual(set(nn[1][0: num[1]]), set([0, 2, 3]))
        self.assertEqual(set(nn[3][0: num[3]]), set([1, 2]))
        nn, num = lc.compute(side, pos, as_array=True, newton=True)
        self.assertEqual(set(nn[0][0: num[0]]), set([1]))
        self.assertEqual(set(nn[1][0: num[1]]), set([0, 2, 3]))
        self.assertEqual(set(nn[3][0: num[3]]), set([1, 2]))
        self.assertTrue(lc.on_border([2.5, 0.0, 0.0]))


if __name__ == '__main__':
    unittest.main()
