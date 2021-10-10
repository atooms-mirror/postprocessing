# Postprocessing

[![pypi](https://img.shields.io/pypi/v/atooms-pp.svg)](https://pypi.python.org/pypi/atooms-pp/)
[![version](https://img.shields.io/pypi/pyversions/atooms-pp.svg)](https://pypi.python.org/pypi/atooms-pp/)
[![license](https://img.shields.io/pypi/l/atooms-pp.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fpostprocessing/HEAD?labpath=docs%2Findex.ipynb)
[![pipeline](https://framagit.org/atooms/postprocessing/badges/master/pipeline.svg)](https://framagit.org/atooms/postprocessing/badges/master/pipeline.svg)
[![coverage report](https://framagit.org/atooms/postprocessing/badges/master/coverage.svg?job=test:f90)](https://framagit.org/atooms/postprocessing/-/commits/master)

Post-processing tools to compute static and dynamic correlation functions from simulations of interacting particles, such as molecular dynamics or Monte Carlo simulations. 

## Quick start

Installation is easy (see [Installation](#installation) for more details)
```
pip install atooms-pp
```

We can now compute correlation functions from trajectories produced
by particle simulation codes. Any trajectory format recognized by
[atooms](https://framagit.org/atooms/atooms.git) can be processed, for instance most "xyz" files
should work fine. If you use a custom trajectory format, it is easy to [add it](https://atooms.frama.io/atooms/).

As an example, we compute the structure factor S(k) for the trajectory
file `trajectory.xyz` contained in the `data/` folder.

### From the command line

![https://www-dft.ts.infn.it/~coslovich/anim.gif](https://framagit.org/atooms/postprocessing/raw/master/docs/anim.gif)

In the example above, we used 20% of the available time frames to compute the averages using the `--norigins` flag. Without it, atooms-pp applies an heuristics to determine the number of time frames required to achieve a reasonable data quality.

The results of the calculation are stored in `data/trajectory.xyz.pp.sk`. If
the system is a mixture of different types of particles, say A and B, the program will create additional files for
partial correlations, named `trajectory.xyz.pp.sk.A-A`, `trajectory.xyz.pp.sk.B-B` and `trajectory.xyz.pp.sk.A-B`.

### From Python

The same calculation can be done from Python:

```python
from atooms.trajectory import Trajectory
import atooms.postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.StructureFactor(t)
     p.do()
```

## Features

Available correlation and distribution functions

- *Real space*
  - radial distribution function
  - mean square displacement
  - velocity auto-correlation function
  - self overlap functions
  - collective overlap functions
  - dynamic susceptibility of the self overlap function
  - non-Gaussian parameter
  - bond-angle distribution
- *Fourier space*
  - structure factor
  - spectral density
  - self intermediate scattering functions
  - collective intermediate scattering functions
  - four-point dynamic susceptibility

## Documentation

Check out the [tutorial](https://atooms.frama.io/postprocessing/index.html) for more examples and the [public API](https://atooms.frama.io/postprocessing/api/postprocessing) for full details.

The tutorial is also available as
- [org-mode file](https://framagit.org/atooms/postprocessing/-/blob/master/docs/index.org)
- [jupyter notebook](https://framagit.org/atooms/postprocessing/-/blob/master/docs/index.ipynb)
- [jupyter notebook on binder](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fpostprocessing/HEAD?labpath=docs%2Findex.ipynb) for interactive execution
- nice, old-school [pdf file](https://framagit.org/atooms/postprocessing/-/blob/master/docs/index.pdf)



## Requirements

- [numpy](https://pypi.org/project/numpy/)
- [atooms](https://framagit.org/atooms/postprocessing.git)
- [optional] [argh](https://pypi.org/project/argh/) (only needed when using `pp.py`)
- [optional] [tqdm](https://pypi.org/project/tqdm/) (enable progress bars)
- [optional] [argcomplete](https://pypi.org/project/argcomplete/) (enable tab-completion for `pp.py`)

## Installation

If you cannot install the package system-wide, you can still install it in the user space. Either from pypi
```
pip install --user atooms-pp
```
or cloning the project repo 
```
git clone https://framagit.org/atooms/postprocessing.git
cd postprocessing
pip install --user .
```
The commands above will install `pp.py` under `~/.local/bin`. Make sure this folder is in your `$PATH`.

## Contributing

Contributions to the project are welcome. If you wish to contribute, check out [these guidelines](https://framagit.org/atooms/atooms/-/blob/master/CONTRIBUTING.md).

## Authors

Daniele Coslovich: http://www-dft.ts.infn.it/~coslovich/
