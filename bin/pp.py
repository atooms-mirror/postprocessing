#!/usr/bin/env python

"""Post processing script."""

import argh
from postprocessing.api import msd, vacf, fkt, fskt, gr, sk, skopti, chi4qs, ik, chi4qs_opti

argh.dispatch_commands([msd, vacf, fkt, fskt, chi4qs, chi4qs_opti, gr, sk, skopti, ik])
