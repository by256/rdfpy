# <img alt="rdfpy" src="./logo.png" height="100">

[![Build Status](https://travis-ci.org/by256/rdfpy.svg?branch=master)](https://travis-ci.org/by256/rdfpy)
[![Documentation Status](https://readthedocs.org/projects/rdfpy/badge/?version=latest)](https://rdfpy.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/rdfpy.svg)](https://pypi.org/project/rdfpy/)
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](https://github.com/by256/rdfpy/blob/master/LICENSE)


**rdfpy** is a Python module for fast computation of 2D and 3D radial distribution functions (RDFs).

## Installation

```shell
$ pip install rdfpy
```

## Usage

```python
import numpy as np
from rdfpy import rdf

# create random particle coordinates in a 20x20x20 box
coords = np.random.uniform(0.0, 20.0, size=(2500, 3))  

# compute radial distribution function with step size = 0.1
g_r, radii = rdf(coords, dr=0.1)
```

You can find a more detailed example in the [Documentation](https://rdfpy.readthedocs.io/).

**Note:** In order for **rdfpy** to work correctly, your particles should spatially be in a cuboidal box, where the entire box is filled with particles.

## How does it work?

**rdfpy** achieves significant speed-up due to:

- **Fast nearest-neighbor look-up**: a k-d tree is utilized when counting the number of particles as a function of distance from an origin particle.
- **Multiprocessing**: computation of the particle count histogram is parallelized across multiple cores, with each core sharing the aforementioned k-d tree.

## Authors

**rdfpy** was developed by [Batuhan Yildirim](https://by256.github.io/) under the supervision of [Prof. Jacqueline M. Cole](https://www.phy.cam.ac.uk/directory/colej).

## Citation

If you use **rdfpy** in your work, please cite:

```
@software{rdfpy,
  author       = {Batuhan Yildirim and
                  Hamish Galloway Brown},
  title        = {by256/rdfpy: rdfpy-v1.0.0},
  month        = mar,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.4625675},
  url          = {https://doi.org/10.5281/zenodo.4625675}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625675.svg)](https://doi.org/10.5281/zenodo.4625675)

## Funding

This project was financially supported by the [Science and Technology Facilities Council (STFC)](https://stfc.ukri.org/) and the [Royal Academy of Engineering](https://www.raeng.org.uk/) (RCSRF1819\7\10).
