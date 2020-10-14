# <img alt="rdfpy" src="./logo.png" height="60">

rdfpy is a Python library for fast computation of 2D and 3D radial distribution functions.

## Installation

```shell
$ pip install rdfpy
```

## Usage

```python
import numpy as np
from rdfpy import rdf3d

particles = np.random.uniform(0.0, 20.0, size=(2500, 3))  # random particles in a 20x20x20 box

g_r, radii = rdf3d(particles, dr=0.1)
```

You can find a more detailed example in the [Documentation](https://rdfpy.readthedocs.io/).

**Note:** In order for rdfpy to work correctly, your particles should spatially be in a cuboidal box, where the entire box is filled with particles.

<!-- 
### What is a Radial Distribution Function?

The radial distribution function (RDF) (or pair correlation function) characterises the structure of a system of particles. If we select an arbitrary particle as the origin, the RDF describes the number of particles we would observe relative to the bulk density of the system, as a function of distance. This is calculated and averaged over every particle in the structure being considered. The formal definition of the RDF is

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=g_{i}(r) = \frac{n_{i}(r)}{4 \pi r^{2}\delta r \rho}">
</p>

where <img src="https://render.githubusercontent.com/render/math?math=n_{i}(r)"> is the number of particles between distances <img src="https://render.githubusercontent.com/render/math?math=r"> and <img src="https://render.githubusercontent.com/render/math?math=r %2B \delta r">, and <img src="https://render.githubusercontent.com/render/math?math=\rho = \frac{N}{V}"> is the number density. Dividing by <img src="https://render.githubusercontent.com/render/math?math=\rho"> ensures that the RDF is centred around 1 when the density of particles observed at some distance does not deviate from the bulk density. Two examples of (1) a highly-ordered crystalline system and (2) a short-range ordered liquid are shown.

<p float="left"  align="center">
  <img src="rdfpy/examples/crystal.png" width="400" />
  <img src="rdfpy/examples/water.png" width="400" /> 
</p> -->

## Authors

[Batuhan Yildirim](https://by256.github.io/)

## Citation

If you use **rdfpy** in your work, please cite:

```
@software{batuhan_yildirim_2020_3932173,
  author = {Batuhan Yildirim},
  title = {by256/rdfpy: rdfpy-v0.1.4},
  month = {July},
  year = {2020},
  publisher = {Zenodo},
  version = {v0.1.4},
  doi = {10.5281/zenodo.3932173},
  url = {https://doi.org/10.5281/zenodo.3932173}
}
```

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3932173.svg)](https://doi.org/10.5281/zenodo.3932173)


