#!/usr/bin/env python

import re
import sys
import argparse
import numpy
from atooms import trajectory
from atooms.utils import add_first_last_skip, fractional_slice
from pyutils.histogram import Histogram


def get_clusters(neighbors, clusters=[]):    
    """ 
    |Neighbor| is a list of lists of neighbors
    |Clusters| is a list of sets of clusters
    """ 
    # TODO: what if we have missing items? Do we expect a list with empty neighbor list?

    my_cluster = {}
    first = 0 

    # Neighbors is empty, return an empty list of clusters
    if len(neighbors) == 0:
        return [[]]

    # No clusters yet, add the first. This iteration could be done before the loop
    if len(clusters) == 0:
        first = 1
        clusters = [set(neighbors[0])]
        for i in neighbors[0]:
            my_cluster[i] = clusters[0]

    # Elements in list are labelled according to the cluster they belong
    # Every time a new element is found is added to label with None as cluster.
    for ne in neighbors[first:]:
        found = None
        # Loop over all elements (neighbors) in ne
        for e in ne:
            for cl in clusters:
                if e in cl:
                    found = cl
                    break
            if found:
                break

        if not found:
            # This is a new cluster
            clusters.append(set(ne))
            for e in ne:
                my_cluster[e] = clusters[-1]
        else:
            # At least one particle belongs to a cluster that is already there
            # We are about to connect to this clusters all elements in ne.
            # We look for each element in ne and if this is already connected to a cluster
            # distinct from the one we just found, we merge the former to the latter.
            for e in ne:
                # Check if this item is connected to previously found clusters
                already_there = e in my_cluster.keys()
                # First loop over all elements in connected clusters and reassign them.
                if already_there:
                    distinct = not my_cluster[e] == found
                    if distinct:
                        # Add former cluster into new one
                        found.update(my_cluster[e])
                        # Remove cluster from list 
                        clusters.remove(my_cluster[e])
                    # Update labels now
                    for ee in my_cluster[e]:
                        my_cluster[ee] = found

                # Now we can update the cluster label of central item
                my_cluster[e] = found

            # Now we add all elements in set
            found.update(ne)

    return clusters

def neighbors_from_voronoi(s, match):
    """Return a list of list of neighbors matching a condition on the central voronoi cell"""
    neigh = []
    for v in s.voronoi:
        if match(v):
            neigh.append([v.central] + list(v.neighbors))
    return neigh

def main(finp, fneigh, ffield, match=None, first=-1, last=-1):

    if fneigh is None:
        fneigh = finp + '.neigh'
    if ffield is None:
        ffield = finp + '.field'
    if first < 0:
        first = 0
    th = trajectory.Sliced(trajectory.Trajectory(finp), slice(first, last, 1))
    tn = trajectory.Sliced(trajectory.TrajectoryNeighbors(fneigh), slice(first, last, 1))
    tf = trajectory.Sliced(trajectory.TrajectoryField(ffield), slice(first, last, 1))

    base = finp
    fout = [base + '.cluster-%s' % match,
            base + '.cluster_largest-%s' % match,
            base + '.cluster_average-%s' % match]

    # Define number of particles in an isolated structure
    hist = Histogram(bin_width=5)

    # Time series and stats of max domain size
    fout = base + '.cluster-%s' % match
    with open(fout, 'w') as fhout:
        fhout.write('# columns: step, average size of %s-clusters, size of largest %s-cluster\n' % 
                    (match, match))
        imin = None
        from atooms import visualize as vis
        from atooms import bonds as bonds
        vis.scene()
        vis.box(th[0].cell)
        vis._radius = {"A": 0.25, "B": 0.15, '0': 0.25, '1': 0.15}
        try:
            match = float(match)
        except:
            pass
        for i in range(len(th)):
            system = th[i]
            step = th.steps[i]
            field = tf[i]
            neighbors_list = []
            for pid, neighbors in enumerate(tn[i].neighbors):
                if imin is None or min(neighbors)<imin:
                    imin =  min(neighbors)
                if field[pid] == match:
                    neighbors_list.append([pid] + 
                                          [ni for ni in neighbors if (field[ni] == match)])
                else:
                    neighbors_list.append([])

            c = get_clusters(neighbors_list)
            x = [len(ci) for ci in c]
            # if len(c) > 0:
            #     hist.add([float(xi) for xi in x])

            step = th.steps[i]
            fhout.write('%d %s %s\n' % (step, numpy.average(x), max(x)))

            # Simply draw all
            #ids = [j for j, p in enumerate(system.particle) if (field[j] == match)]
            #vis.show([system.particle[j] for j in ids], radius_auto=False)
            cluster_particle = []
            for ids in c:
                if len(ids) == 0:
                    continue
                root = system.particle[list(ids)[0]]
                cluster_particle += [system.particle[j].nearest_image(root, system.cell) for j in ids]
                bonds.subset_with_bonds(ids, system.particle, cell=system.cell, neighbors=neighbors_list)
            vis.show(cluster_particle, radius_auto=False)
            vis.pause()
            bonds.clear()

        fhout.write(hist.stats)

    # Histogram of domain size
    # with open(base + '.cluster-%s.stats' % match, 'w') as fh:
    #     fh.write(hist.stats)
    #     fh.write(str(hist))

if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
                           
#     # parser = OptionParser(usage='%prog [options] <args>')
#     # parser.add_option('--verbose', action='store_true', dest='', help='verbose mode')
#     # opts, args = parser.parse_args()
#     ClusterAnalysis([[1,2,3], [4,5,6], [5,8,9], [10,11], [11,12], [1,12]])

