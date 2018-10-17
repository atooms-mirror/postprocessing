Post processing
==================
Python post processing tools to compute static and dynamic correlation functions from particle simulations

- Real space: radial distribution function, mean square displacement, time-dependent overlap functions, non-Gaussian parameter
- Fourier space: structure factor, intermediate scattering functions, four-point dynamic susceptibility

This package relies on [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git) to read trajectory files.

Quick start
------------
Installation is easy (see [Installation] for more details)
```
pip install atooms-pp
```

We can now compute correlation functions from trajectories produced
by particle simulation codes. Any trajectory format recognized by
atooms can be processed, for instance most "xyz" files
should work fine.

As an example, we compute the structure factor S(k) for the trajectory
file `trajectory.xyz` contained in the `data/` directory.

![https://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/anim.gif](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing/raw/develop/docs/anim.gif)

In the example above, we used 20% of the available time frames to compute the averages using the `--norigins` flag. Without it, atooms-pp applies an heuristics to determine the number of time frames required to achieve a reasonable data quality.

The results of the calculation are stored in `data/trajectory.xyz.pp.sk`. If
the system is a mixture of different types of particles, say A and B, the program will create additional files for
partial correlations, named `trajectory.xyz.pp.sk.A-A`, `trajectory.xyz.pp.sk.B-B` and `trajectory.xyz.pp.sk.A-B`.

The same calculation can be done from python:

```python
from atooms.trajectory import Trajectory
import atooms.postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.StructureFactor(t)
     p.do()
```

Checkout the [tutorial](https://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/pp_notebook/) and [notebook](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing/raw/develop/docs/tutorial.ipynb) for more details.

Requirements
------------
- [numpy](https://pypi.org/project/numpy/)
- [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git)
- [optional] [argh](https://pypi.org/project/argh/) (only needed when using `pp.py`)
- [optional] [tqdm](https://pypi.org/project/tqdm/) (enable progress bars)
- [optional] [argcomplete](https://pypi.org/project/argcomplete/) (enable tab-completion for `pp.py`)

Installation
------------
If you cannot install the package system-wide, you can still install it in the user space. Either from pypi
```
pip install --user atooms-pp
```
or cloning the project repo 
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git
cd postprocessing
make user
```
The commands above will install `pp.py` under `~/.local/bin`. Make sure this folder is in your `$PATH`. To install system-wide, `sudo make install`.