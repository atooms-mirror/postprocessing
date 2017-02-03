Post processing
==================

Python analysis tools for particle-based simulations.

- Static and time-dependent correlation functions
  - real space: radial distribution function, mean square displacement, time-dependent overlap, ...
  - Fourier space: structure factor, intermediate scattering function, dynamic susceptibility, ...
- Bond-orientational order
- Common neighbors analysis

This package relies on [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git) to read trajectory files.

Getting started
---------------

Post processing tools typically operate on trajectory files produced
by molecular simulation codes. Any trajectory format recognized by
atooms can be processed, for instance most "xyz" trajectory files should work fine. 
Most tools are simple python scripts that can
be executed from the command line. 

Example: the following command
will compute the radial distribution function g(r) from the trajectory
file `trajectory.xyz` contained in the `data/` directory

```bash
$ pp.py gr data/trajectory.xyz
```

The results will be stored in the file `data/trajectory.xyz.pp.gr`. If
multiple chemical species are present, the program will create files
named `trajectory.xyz.pp.gr.1-1`, `trajectory.xyz.pp.gr.2-2` and so on.

The same kind of calculation can be done from python:

```python
from atooms.trajectory import Trajectory
import postprocessing as pp

with Trajectory('data/trajectory.xyz') as t:
     p = pp.RadialDistributionFunction(t)
     p.do()
```

Requirements
------------
- numpy
- [atooms](https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git)

Installation
------------
From the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing.git
cd postprocessing
make install
```