#!/usr/bin/env python

"""Post processing script."""

import argh
from postprocessing.api import msd, vacf, fkt, fskt, gr, sk

argh.dispatch_commands([msd, vacf, fkt, fskt, gr, sk])
