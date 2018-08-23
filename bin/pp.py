#!/usr/bin/env python

"""Post processing script."""

import argh
import argparse
from atooms.core.utils import setup_logging
from atooms.core.utils import add_first_last_skip
import atooms.postprocessing as postprocessing
from atooms.postprocessing.api import msd, vacf, fkt, fskt, gr, sk, skopti, chi4qs, ik, chi4qs_opti


# We add some global some flags. For backward compatibility, we keep
# them in the function signature as well.
parser = argparse.ArgumentParser()
parser = add_first_last_skip(parser)
parser.add_argument('--fmt', dest='fmt', help='fmt')
parser.add_argument('--output', dest='output', default='{trajectory.filename}.pp.{short_name}.{tag}', help='output path like {pp.trajectory.filename}.pp.{pp.short_name}.{pp.tag}')
parser.add_argument('--verbose', action='store_true', dest='verbose', help='verbose output')
parser.add_argument('--debug', action='store_true', dest='debug', help='debug output')
parser.add_argument('--species-layout', dest='species_layout', help='force species layout to F, C or A')
argh.add_commands(parser, [msd, vacf, fkt, fskt, chi4qs, chi4qs_opti, gr, sk, skopti, ik])
args = parser.parse_args()

postprocessing.correlation.OUTPUT_PATH = args.output 
if args.verbose:
    setup_logging('postprocessing', level=20)
    setup_logging('atooms', level=40)
elif args.debug:
    setup_logging('postprocessing', level=10)
    setup_logging('atooms', level=40)
else:
    setup_logging('postprocessing', level=40)
    setup_logging('atooms', level=40)

argh.dispatching.dispatch(parser)
