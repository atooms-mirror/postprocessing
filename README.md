Post processing
==================

Tools to analyze molecular simulation data. Current

- Time dependent correlation functions
  - real space: radial distribution function, mean square displacement, time-dependent overlap 
  - Fourier space: structure factor, intermediate scattering function, dynamic susceptibility
- Bond-orientational order
- Common neighbors analysis

Getting started
---------------

Post-processing tools typically take as input trajectory files
produced by some molecular simulation code. This can be done from the
command line. For instance, the following line will compute the radial
distribution function g(r) from the trajectory file `trajectory.xyz`

```bash
$ pp.py gr trajectory.xyz
```

Any trajectory format supported by `atooms` can be passed to the post-processing tools.

The same calculation can be done from python

```python
from atooms.trajectory import Trajectory
import postprocessing as pp

with Trajectory('trajectory.xyz') as t:
     p = pp.RadialDistributionFunction(t)
     p.do()
```

Requirements
------------
- numpy
- atooms

Installation
------------
From the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocession.git
make install
```
