# rdfpy

rdfpy is a Python package for computing 2D and 3D radial distribution functions. 

![Alt text](rdfpy/examples/random2d.png)

## Installation

```shell
$ pip install rdfpy
```

## Usage

```python
import numpy as np
from rdfpy import rdf3d

r = 0.5  # arbitrary radius for random example
particles = np.random.uniform(0.0, 10.0, size=(1000, 3))  # random particles in a 10x10x10 box

g_r, radii = rdf3d(particles, r, h=10, w=10, d=10, n_radii=128)
```

## Authors

[Batuhan Yildirim](http://www.mole.phy.cam.ac.uk/people/by.php)

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)
