#!/usr/bin/env python
import sys
import atooms.trajectory as trj
import atooms.voronoi as vor

def bond(f, broken_bonds=1):
#    with trj.TrajectoryXYZ(f) as th:
    with vor.TrajectoryVoronoi(f) as th:
        npart = len(th[0].particle)
        
        # Put data into a big list
        from collections import namedtuple
        voronoi_list = []
        for j in range(npart):
            voronoi_list.append([])
        for frame, s in enumerate(th):
            for i in range(npart):
                #nn = sorted([int(_) for _ in s.particle[i].neighbors.split(',')])
                #nn = sorted(s.particle[i].neighbors)
                nn = set(s.particle[i].neighbors)
                #nn = sorted([int(_) for _ in getattr(s.particle[i], 'neighbors*').split(',')])
                voronoi_list[i].append(nn)
                
                #        for t in voronoi_list[820]:
                #            print t
        
        from collections import defaultdict
        prob = defaultdict(list)
        for i in range(npart):
            for t in [0, 1, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 92, 128]:
                for t0 in range(0, len(voronoi_list[i])-t, 10):
                    #                    print t, t0, voronoi_list[i][t0], voronoi_list[i][t0+t], voronoi_list[i][t0] == voronoi_list[i][t0+t]
                    # prob[t].append(int(voronoi_list[i][t0] == voronoi_list[i][t0+t]))
                    prob[t].append(len(voronoi_list[i][t0] & voronoi_list[i][t0+t]) > broken_bonds)

        for key in sorted(prob):
            print key, sum(prob[key]) / float(len(prob[key])) #, len(prob[key])
#bond('/home/coslo/usr/atooms/voronoi/data/hspoly.h5.voronoi.xyz')

import argh
argh.dispatch_command(bond)

