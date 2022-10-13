#!/usr/bin/env python

import sys
import os
import random
import unittest
import numpy
from atooms.postprocessing import helpers


def deviation(x, y):
    return (numpy.sum((x-y)**2)/len(x))**0.5

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_ifabsmm(self):
        x = numpy.array(numpy.linspace(0.0, 1.0, 20))
        f = (x-0.75)**2
        results = helpers.ifabsmm(x, f)
        xmin, fmin = results[0]
        xmax, fmax = results[1]
        self.assertAlmostEqual(xmin, 0.75)
        self.assertAlmostEqual(xmax, 0.0)

        f = (x-1.0)**2
        results = helpers.ifabsmm(x, f)
        xmin, fmin = results[0]
        self.assertAlmostEqual(xmin, 1.0)


if __name__ == '__main__':
    unittest.main()
