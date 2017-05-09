#!/usr/bin/env python

"""Bond-orientational order parameters."""

import os
import sys
import numpy
from atooms.trajectory import Trajectory, TrajectoryXYZ
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
    """
    Average q_l as a function of time. Values of l must be passed as
    args.lvalues, a string containing comma separated entries,
    Example: lvalues='4,6'.
    """

    t = Trajectory(f)
    tn = get_neighbors(f, None, args) #, os.path.basename(sys.argv[0]))

    if args.field_file is not None:
        tf = TrajectoryXYZ(args.field_file)
        if args.field is None:
            raise ValueError('provide field')

    desc = ''
    mode = 'w'
    fbase = f + '.boo%s' % args.tag
    lvalues = [int(i) for i in args.lvalues.split(',')]

    fq, fqb = {}, {}
    for l in lvalues:
        # Keep track of which neighbor file was used
        fq[l] = open(fbase + '.q%d' % l, 'w', buffering=0)
        fq[l].write('# neighbors: %s\n' % desc)
        # If requested, write Lechner-Dellago variants too
        if not args.nobar:
            fqb[l] = open(fbase + '.qb%d' % l, 'w', buffering=0)
            fqb[l].write('# neighbors: %s\n' % desc)

    for i, s in enumerate(t):
        step = t.steps[i]
        # Find sample that matches step in neighbor file
        # If not found, skip sample.
        try:
            index_neigh = tn.steps.index(step)
        except:
            continue

        if args.field_file is not None:
            index_field = tf.steps.index(step)
            field = []
            for pi in tf[index_field].particle:
                fi = getattr(pi, args.field)
                fi = [float(x) for x in fi.split(',')]
                field.append(fi)
        else:
            field = None

        # Compute average bond orientational order
        print 'boo compute step', step, '...',
        sys.stdout.flush()
        b = boo.BondOrientationalOrder(s.particle, tn[index_neigh].neighbors, s.cell.side, field)
        q, qb = {}, {}
        for l in lvalues:
            q[l] = b.ql(l)
            if not args.nobar:
                qb[l] = b.ql_bar(l)
        print 'done'

        # Write in xyz format
        if args.xyz:
            tag = 'q'.join([str(i) for i in lvalues])
            cols = ',q'.join([str(i) for i in lvalues])
            cols = cols[1:]
            # TODO: we can use trajectory field here
            write_xyz(fbase + '.%s.xyz' % tag, [q[l] for l in lvalues], {'step':step, 'columns':cols}, mode)
            mode = 'a'
            
        # Dump average
        for l in lvalues:
            fq[l].write('%d %g\n' % (step, numpy.average(q[l])))
            if not args.nobar:
                fqb[l].write('%d %g\n' % (step, numpy.average(qb[l])))

    for l in lvalues:
        fq[l].close()
        if not args.nobar:
            fqb[l].close()

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
    parser.add_argument(      '--field',   dest='field', default='field', help='field to read')
    parser.add_argument('-F', '--field-file', dest='field_file', default=None, help='field file in xyz format')
    parser.add_argument('-l', '--l-values', dest='lvalues', default='4,6', help='comma separated list of l (el) values')
    parser.add_argument(nargs='+',         dest='files',type=str, help='input files')
    args = parser.parse_args()

    if len(args.tag) > 0:
        args.tag = '-' + args.tag

    for f in args.files:
        if args.cluster:
            ave_cluster(f, args)
        else:
            if args.map:
                map_boo(f, args)
            else:
                ave(f, args)
