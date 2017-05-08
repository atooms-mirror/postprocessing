#!/usr/bin/env python

"""Bond-orientational order parameters."""

import os
import sys
import numpy
from atooms.trajectory import Trajectory
import postprocessing.boo as boo
from postprocessing.neighbors import all_neighbors, get_neighbors

def write_xyz(filename, data, meta, mode='w'):
    # TODO: move it up the chain
    # Check that all arrays in data have the same length
    nlines = len(data[0])
    ncols = len(data)
    lengths_ok = map(lambda x: len(x) == nlines, data)
    if not all(lengths_ok):
        raise ValueError('All arrays must have the same length')

    # Write in xyz format    
    with open(filename, mode) as fh:
        fh.write('%d\n' % nlines)
        # Comment line: concatenate metadata
        line = ''
        for key in meta:
            line += '%s:%s ' % (key, meta[key])
        line += '\n'
        fh.write(line)
        # Write data. This is inefficient but
        # we cannot use numpy.savetxt because there is no append mode.
        fmt = '%s ' * len(data)
        fmt = fmt[:-1] + '\n'
        for i in range(nlines):
            fh.write(fmt % tuple([data[j][i] for j in range(ncols)]))

def ave(f, args):
    """Average q4 and q6 as a function of time.
    Dump map if requested."""

    t = Trajectory(f)
    tn = get_neighbors(f, None, args) #, os.path.basename(sys.argv[0]))
    if args.field_file is not None:
        tf = TrajectoryField(args.field_file)
        if args.field is None:
            raise ValueError('provide field')

    desc = ''
    mode = 'w'
    fbase = f + '.boo%s' % args.tag

    fq4 = open(fbase + '.q4', 'w', buffering=0)
    fq6 = open(fbase + '.q6', 'w', buffering=0)
    # Keep track of which neighbor file was used
    fq4.write('# neighbors: %s\n' % desc)
    fq6.write('# neighbors: %s\n' % desc)

    # If requested, write Lechner-Dellago variants too
    if not args.nobar:
        fqb4 = open(fbase + '.qb4', 'w', buffering=0)
        fqb6 = open(fbase + '.qb6', 'w', buffering=0)
        # Keep track of which neighbor file was used
        fqb4.write('# neighbors: %s\n' % desc)
        fqb6.write('# neighbors: %s\n' % desc)

    for i, s in enumerate(t):
        step = t.steps[i]
        # Find sample that matches step in neighbor file
        # If not found, skip sample.
        try:
            index_neigh = tn.steps.index(step)
        except:
            continue
        if args.field is not None:
            index_field = tf.steps.index(step)
            field = getattr(tf[index_field].particle, args.field)
        else:
            field = None

        # Compute average bond orientational order
        print 'boo compute step', step, '...',
        sys.stdout.flush()
        field = None
        b = boo.BondOrientationalOrder(s.particle, tn[j].neighbors, s.cell.side, field)
        q4 = b.ql(4)
        q6 = b.ql(6)
        if not args.nobar:
            qb4 = b.ql_bar(4)
            qb6 = b.ql_bar(6)
        print 'done'

        # Write q4, q6 in xyz format
        if args.xyz:
            write_xyz(fbase + '.q4q6.xyz', [q4, q6], {'step':step, 'columns':'q4,q6'}, mode)
            mode = 'a'
            
        # Dump average
        fq4.write('%d %g\n' % (step, numpy.average(q4)))
        fq6.write('%d %g\n' % (step, numpy.average(q6)))
        if not args.nobar:
            fqb4.write('%d %g\n' % (step, numpy.average(qb4)))
            fqb6.write('%d %g\n' % (step, numpy.average(qb6)))

    fq4.close()
    fq6.close()
    if not args.nobar:
        fqb4.close()
        fqb6.close()

    if len(args.neigh_file) == 0:
        os.remove(tn.filename)


def map_boo(f, args):
    """Map q4 and q6 as a function of time."""

    t = Trajectory(f)
    tn = get_neighbors(f, None, args) #get_neighbors(f, args, os.path.basename(sys.argv[0]))

    # This is a full dump on a per-particle basis
    if not args.nobar:
        fname = f + '.boo%s.qb4qb6' % args.tag
    else:
        fname = f + '.boo%s.q4q6' % args.tag

    with open(fname, 'w') as fmap:
        for i, s in enumerate(t):
            if i % args.skip != 0:
                continue
            step = t.steps[i]
            # Find sample that matches step in neighbor file
            # If not found, skip sample.
            try:
                j = tn.steps.index(step)
            except:
                continue
            # Compute average bond orientational order
            b = boo.BondOrientationalOrder(s.particle, tn[j].neighbors, s.cell.side)
            if not args.nobar:
                q4 = b.ql_bar(4)
                q6 = b.ql_bar(6)
            else:
                q4 = b.ql(4)
                q6 = b.ql(6)

            # Dump q4,q6 map
            fmap.write('#\n' % step)
            fmap.write('# step:%d neighbors: %s\n' % (step, ''))
            for i, j in zip(q4, q6):
                fmap.write('%g %g\n' % (i, j))

    if len(args.neigh_file) == 0:
        os.remove(tn.filename)


def ave_cluster(f, args):

    """BOO for a cluster."""

    t = Trajectory(f)

    with open(f + '.boo%s.q4' % args.tag, 'w', buffering=0) as fq4,\
         open(f + '.boo%s.q6' % args.tag, 'w', buffering=0) as fq6:

        for i, s in enumerate(t):
            step = t.steps[i]
            # Compute average bond orientational order
            b = boo.BondOrientationalOrder(s.particle)
            q4 = b.ql(4)
            q6 = b.ql(6)
            fq4.write('%d %g\n' % (step, q4[0]))
            fq6.write('%d %g\n' % (step, q6[0]))

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-n',              dest='neigh_file', type=str, default='', help='neighbors file')
    parser.add_argument('-N', '--neighbor',dest='neigh', type=str, default='', help='flags for neigh.x command')
    parser.add_argument('-V', '--neighbor-voronoi',dest='neigh_voronoi', action='store_true', help='neigh_file is of Voronoi type')
    parser.add_argument('-M', '--neighbor-max',dest='neigh_limit', type=int, default=None, help='take up to *limit* neighbors (assuming they are ordered)')
    parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
    parser.add_argument('-m', '--map',     dest='map', action='store_true', help='only q4,q6 map (no averages)')
    parser.add_argument(      '--xyz',     dest='xyz', action='store_true', help='write q4,q6 in xyz format')
    parser.add_argument('-B', '--no-bar',  dest='nobar', action='store_true', help='do not compute Lechner-Dellago variants')
    parser.add_argument(      '--cluster', dest='cluster', action='store_true', help='only compute boo for first particles')
    parser.add_argument(      '--skip',    dest='skip', default=1, help='skip every SKIP configuration')
    parser.add_argument(nargs='+',         dest='files',type=str, help='input files')
    args = parser.parse_args()

    if len(args.tag) > 0:
        args.tag = '_' + args.tag

    for f in args.files:
        if args.cluster:
            ave_cluster(f, args)
        else:
            if args.map:
                map_boo(f, args)
            else:
                ave(f, args)
