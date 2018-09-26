Post processing
==================

Python post processing tools to compute static and dynamic correlation functions from particle-based simulations. Supported correlation functions:
- real space: radial distribution function, mean square displacement, time-dependent overlap, ...
- Fourier space: structure factor, intermediate scattering function, dynamic susceptibility, ...

This package relies on [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git) to read trajectory files.

Getting started
---------------

Post processing tools typically operate on trajectory files produced
by molecular simulation codes. Any trajectory format recognized by
atooms can be processed, for instance most "xyz" trajectory files
should work fine. The correlation functions can be computed using
either the command line script `pp.py` or directly from python.

Example: the following command
will compute the radial distribution function g(r) from the trajectory
file `trajectory.xyz` contained in the `data/` directory

```bash
$ pp.py gr data/trajectory.xyz
```

The results will be stored in the file `data/trajectory.xyz.pp.gr`. If
multiple chemical species are present, say A and B, the program will create additional files for
partial correlations, named `trajectory.xyz.pp.gr.A-A`, `trajectory.xyz.pp.gr.B-B` and `trajectory.xyz.pp.gr.A-B`.

The same kind of calculation can be done from python:

```python
from atooms.trajectory import Trajectory
import postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.RadialDistributionFunction(t)
     p.do()
```

See the tutorials under `docs/` for more details.


Requirements
------------
- numpy
- [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git)
- argh (optional, only needed when using `pp.py`)

Installation
------------
From the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git
cd postprocessing
make install
```
