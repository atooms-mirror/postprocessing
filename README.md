Post processing
==================

Python post processing tools to compute static and dynamic correlation functions from particle simulations. Supported correlation functions:
- real space: radial distribution function, mean square displacement, time-dependent overlap, ...
- Fourier space: structure factor, intermediate scattering function, dynamic susceptibility, ...

This package relies on [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git) to read trajectory files.

Getting started
---------------

Our goal is to analyze trajectories produced
by molecular simulation codes. Any trajectory file format recognized by
atooms can be processed, for instance most "xyz" trajectory files
should work fine. The correlation functions can be computed using
either the command line script `pp.py` or directly from python.

![terminal](docs/anim.gif)

Example: we compute the radial distribution function g(r) for the trajectory
file `trajectory.xyz` contained in the `data/` directory

```bash
$ pp.py gr data/trajectory.xyz
```

The results are stored in `data/trajectory.xyz.pp.gr`. If
multiple particle species are present, say A and B, the program will create additional files for
partial correlations, named `trajectory.xyz.pp.gr.A-A`, `trajectory.xyz.pp.gr.B-B` and `trajectory.xyz.pp.gr.A-B`.

The same kind of calculation can be done from python:

```python
from atooms.trajectory import Trajectory
import atooms.postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.RadialDistributionFunction(t)
     p.do()
```

See the tutorials under `docs/` for more details.

Requirements
------------
- numpy
- [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git)
- [optional] argh (only needed when using `pp.py`)
- [optional] tqdm (only needed to enable progress bars)

Installation
------------
From the python package index
```
pip install atooms-pp
```

From the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git
cd postprocessing
make user
```
