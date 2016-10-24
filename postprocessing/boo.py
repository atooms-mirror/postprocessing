# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import time
import numpy
import copy
from scipy.special import sph_harm as Y

def periodic_vector(a, box):
    for i in xrange(len(box)):
        a[...,i] = numpy.where(a[...,i] >  box[i]/2, a[...,i] - box[i], a[...,i])
        a[...,i] = numpy.where(a[...,i] < -box[i]/2, a[...,i] + box[i], a[...,i])
    return a

def car2sph(rvec):
    # from http://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
    sp = numpy.ndarray(rvec.shape)
    xy = rvec[:,0]**2 + rvec[:,1]**2
    sp[:,0] = numpy.sqrt(xy + rvec[:,2]**2)
    sp[:,1] = numpy.arctan2(rvec[:,1], rvec[:,0]) # longitude
    sp[:,2] = numpy.arctan2(numpy.sqrt(xy), rvec[:,2]) # latitude
    return sp

import warnings
warnings.simplefilter("ignore", numpy.ComplexWarning)

class BondOrientationalOrder(object):
    
    def __init__(self, particle, neighbors=None, box=None):
        """If neighbors and box are None, this is a cluster and we compute
        a full neighbor list around the central particle."""
        self.position = numpy.array([p.position for p in particle])
        self.neighbors = neighbors
        self.box = box

    def _qlm(self, l):
        # TODO: optimize kernel by packing more particles in multi dim array and passing them to Y. Y is waht takes time
        np = self.position.shape[0]
        qlm = numpy.zeros((2*l+1, np), dtype=complex)
        for i in xrange(np):            
            #nn = numpy.asarray(self.neighbors[i][:])
            nn = numpy.asarray(self.neighbors[i][:])
            rvec = self.position[i, :] - self.position[nn, :]
            rvec = periodic_vector(rvec, self.box)
            sph = car2sph(rvec)
            for m in xrange(2*l + 1):
                qlm[m, i] = numpy.average(Y(m-l, l, sph[:,1], sph[:,2]))
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
            qlm[m, i] = numpy.average(Y(m-l, l, sph[:,1], sph[:,2]))
        return qlm

    def _rot_invariant(self, q, l):
        """ Construct rotational invariant or order l from full boo parameter """
        s = numpy.sum(q*q.conj(), axis=0)
        s = numpy.array(s, dtype=float)
        return numpy.sqrt(4*numpy.pi / (2*l+1) * s)

    def ql(self, l=6):
        # Check whether this is a cluster
        if self.neighbors is None:
            return self._rot_invariant(self._qlm_cluster(l), l)
        else:
            return self._rot_invariant(self._qlm(l), l)

    def ql_bar(self, l=6):
        qlm = self._qlm(l)
        qlm_bar = numpy.zeros_like(qlm)        
        for i in xrange(qlm.shape[1]):
            nn = self.neighbors[i][:]
            qlm_bar[:, i] = qlm[:, i] + numpy.sum(qlm[:, nn], axis=1)
            #qlm_bar[:, i] = qlm[:, i] + numpy.sum(qlm[:, nn[0:min(len(nn), 12)]], axis=1)
            qlm_bar[:, i] /= len(self.neighbors[i]) + 1
            #qlm_bar[:, i] /= min(len(nn), 12) + 1
        return self._rot_invariant(qlm_bar, l)

#TODO: refactor sample looping using callbacks
def dump(trajectory, neighbors, fmt='xyz'):

    s = trajectory.read(0)
    box = s.cell.side
    np = len(s.particle)

    import atooms.trajectory as trj

    if fmt[-3:]=='xyz':
        t = trj.TrajectoryXYZ(trajectory.filename + '.boo.' + fmt, 'w')
    else:
        t = trj.TrajectoryPDB(trajectory.filename + '.boo.pdb', 'w')
    t.write_initial_state(s)

    for i, step in zip(trajectory.samples, trajectory.steps):
        try:
            j = neighbors.steps.index(step)
        except:
            print 'could not find %d in neighbors' % step
            continue
        neighbors.samples[j]
        s = trajectory.read(i)
        n = neighbors.read(j)
        boo = BondOrientationalOrder(s.particle, n, box)
        for p, qi in zip(s.particle, boo.ql_bar(6)):
            p.tag = qi
        t.write_sample(s, sample=i, step=step)
 
    t.close()
        
def analyze(trajectory, neighbors, nbins=20):

    from pyutils.histogram import histogram, histogram2d

    s = trajectory.read(0)
    box = s.cell.side
    np = len(s.particle)

    q6 = []
    q4 = []
    for i, step in zip(trajectory.samples, trajectory.steps):
        try:
            j = neighbors.steps.index(step)
        except:
            print 'could not find %d in neighbors' % step
            continue
        neighbors.samples[j]
        s = trajectory.read(i)
        n = neighbors.read(j)
        boo = BondOrientationalOrder(s.particle, n, box)
        q4 += list(boo.ql_bar(4))
        q6 += list(boo.ql_bar(6))

    histograms = {'q4' : histogram(q4, nbins),
                  'q6' : histogram(q6, nbins), 
                  'q46': histogram2d(q4, q6, nbins)
                  }

    averages = {'q4' : numpy.average(q4),
                'q6' : numpy.average(q6),
                'q46' : 0.0
                }

    return averages, histograms

def write(fname, ave, hist):
    import pyutils.table
    for q in ['q4', 'q6', 'q46']:
        with open('%s.boo.%s' % (fname, q), 'w') as fh:
            fh.write('# ave = %g\n' % ave[q])
            fh.write(pyutils.table.table(hist[q]))

#from atooms import trajectory as trj
#d = '/media/data/projects/cavity_hocky/data/wca-poly_P_N1000_T0.03_n00/p0.58000/config.dat.lin'
#h=dump(trj.TrajectoryHDF5(d), trj.TrajectoryNeighbors(d+'.voronoi.xyz.neigh'))
#print analyze(trj.TrajectoryHDF5(d), trj.TrajectoryNeighbors(d+'.voronoi.xyz.neigh'))

#def plot():
    # pylab.imshow(h, extent=(x[0], x[-1], y[0], y[-1]), interpolation='nearest')
    # pylab.colorbar()
    # pylab.show()
