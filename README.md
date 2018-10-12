Post processing
==================

Python post processing tools to compute static and dynamic correlation functions from particle simulations. Supported correlation functions:
- real space: radial distribution function, mean square displacement, time-dependent overlap, ...
- Fourier space: structure factor, intermediate scattering function, dynamic susceptibility, ...

This package requires [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git) to read trajectory files. It supports tab-completion of command line arguments and progress bar via the optional packages `argcomplete` and `tqdm`.

Quick start
---------------
Our goal is to analyze trajectories produced
by particle simulation codes. Any trajectory file format recognized by
atooms can be processed, for instance most "xyz" trajectory files
should work fine. The correlation functions can be computed using
either the command line script `pp.py` or directly from python.

Here we compute the structure factor S(k) for the trajectory
file `trajectory.xyz` contained in the `data/` directory.

![terminal](docs/anim.gif)

In the example above, we actually used 20% of the available
time frames ("time origins"). If this is not specified, atooms-pp applies an heuristics to determine the number of time frames needed to achieve a reasonable data quality.

```bash
$ pp.py sk data/trajectory.xyz
```

The results of the calculations are stored in `data/trajectory.xyz.pp.sk`. If
multiple particle species are present, say A and B, the program will create additional files for
partial correlations, named `trajectory.xyz.pp.sk.A-A`, `trajectory.xyz.pp.sk.B-B` and `trajectory.xyz.pp.sk.A-B`.

The same calculation can be done from python:

```python
from atooms.trajectory import Trajectory
import atooms.postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.StructureFactor(t)
     p.do()
```

See the tutorials under `docs/` for more details.

Requirements
------------
- numpy
- [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git)
- [optional] argh (only needed when using `pp.py`)
- [optional] tqdm (only needed to enable progress bars)
- [optional] argcomplete (only needed to enable tab-completion for `pp.py`)

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
This will install `pp.py` under `~/.local/bin`. Make sure this folder is in your `$PATH`. Alternatively, install the package system-wide with `sudo make install`.
