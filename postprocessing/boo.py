# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Bond-orientational order parameters."""

import numpy
from scipy.special import sph_harm
import warnings
warnings.simplefilter("ignore", numpy.ComplexWarning)

# Helper functions

def car2sph(rvec):
    """Cartesian to spherical coordinates conversion.

    Source: http://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
    """
    sp = numpy.ndarray(rvec.shape)
    xy = rvec[:,0]**2 + rvec[:,1]**2
    sp[:,0] = numpy.sqrt(xy + rvec[:,2]**2)
    sp[:,1] = numpy.arctan2(rvec[:,1], rvec[:,0]) # longitude
    sp[:,2] = numpy.arctan2(numpy.sqrt(xy), rvec[:,2]) # latitude
    return sp

def periodic_vector(a, box):
    for i in xrange(len(box)):
        a[...,i] = numpy.where(a[...,i] >  box[i]/2, a[...,i] - box[i], a[...,i])
        a[...,i] = numpy.where(a[...,i] < -box[i]/2, a[...,i] + box[i], a[...,i])
    return a


class BondOrientationalOrder(object):
    
    def __init__(self, particle, neighbors=None, box=None):
        """If `neighbors` and `box` are None, particles are treated as an
        isolated cluster and we consider all particles around the
        central one as neighbors.
        """
        self.position = numpy.array([p.position for p in particle])
        self.neighbors = neighbors
        self.box = box

    def _qlm(self, l):
        np = self.position.shape[0]
        qlm = numpy.zeros((2*l+1, np), dtype=complex)
        for i in xrange(np):            
            nn = numpy.asarray(self.neighbors[i][:])
            rvec = self.position[i, :] - self.position[nn, :]
            rvec = periodic_vector(rvec, self.box)
            sph = car2sph(rvec)
            for m in xrange(2*l + 1):
                Y = sph_harm(m-l, l, sph[:,1], sph[:,2])
                qlm[m, i] = numpy.average(Y)
        return qlm

    def _qlm_cluster(self, l):
        # TODO: get particle closest to CM
        np = self.position.shape[0]
        qlm = numpy.zeros((2*l+1, 1), dtype=complex)
        i = 0
        # In a cluster, we get all neighbors 
        nn = numpy.array([j for j in range(np) if i != j])
        rvec = self.position[i, :] - self.position[nn, :]
        sph = car2sph(rvec)
        for m in xrange(2*l + 1):
            qlm[m, i] = numpy.average(sph_harm(m-l, l, sph[:,1], sph[:,2]))
        return qlm

    def _rot_invariant(self, q, l):
        """Construct rotational invariant or order l from full boo parameter."""
        s = numpy.sum(q*q.conj(), axis=0)
        s = numpy.array(s, dtype=float)
        return numpy.sqrt(4*numpy.pi / (2*l+1) * s)

    def ql(self, l=6):
        if self.neighbors is None:
            # This is an isolated cluster
            return self._rot_invariant(self._qlm_cluster(l), l)
        else:
            return self._rot_invariant(self._qlm(l), l)

    def ql_bar(self, l=6):
        """Lechner-Dellago variant."""
        qlm = self._qlm(l)
        qlm_bar = numpy.zeros_like(qlm)        
        for i in xrange(qlm.shape[1]):
            nn = self.neighbors[i][:]
            qlm_bar[:, i] = qlm[:, i] + numpy.sum(qlm[:, nn], axis=1)
            qlm_bar[:, i] /= len(self.neighbors[i]) + 1
            #qlm_bar[:, i] = qlm[:, i] + numpy.sum(qlm[:, nn[0:min(len(nn), 12)]], axis=1)
            #qlm_bar[:, i] /= min(len(nn), 12) + 1
        return self._rot_invariant(qlm_bar, l)

