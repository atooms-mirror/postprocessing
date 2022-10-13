#!/usr/bin/env python

import sys
import os
import random
import unittest
import numpy
import atooms.postprocessing.api as api


def deviation(x, y):
    return (numpy.sum((x-y)**2)/len(x))**0.5

class TestApi(unittest.TestCase):

    def setUp(self):
        self.reference_path = 'data'
        if not os.path.exists(self.reference_path):
            self.reference_path = os.path.join(os.path.dirname(sys.argv[0]), '../data')
        self.test_file = os.path.join(self.reference_path, 'kalj-small.xyz')

    def test_gr(self):
        api.gr(self.test_file)

    def test_sk(self):
        api.sk(self.test_file)

    def test_ik(self):
        api.ik(self.test_file)

    def test_msd(self):
        api.msd(self.test_file)

    def test_vacf(self):
        api.vacf(self.test_file)

    def test_fskt(self):
        api.fskt(self.test_file, dk=0.2)

    def test_fkt(self):
        api.fkt(self.test_file, dk=0.2)

    def test_chi4qs(self):
        api.chi4qs(self.test_file)


if __name__ == '__main__':
    unittest.main()
