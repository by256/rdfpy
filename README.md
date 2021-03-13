# <img alt="rdfpy" src="./logo.png" height="60">

**rdfpy** is a Python library for fast computation of 2D and 3D radial distribution functions.

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

**Note:** In order for **rdfpy** to work correctly, your particles should spatially be in a cuboidal box, where the entire box is filled with particles.

## Authors

**rdfpy** was developed by [Batuhan Yildirim](https://by256.github.io/) under the supervision of [Prof. Jacqueline M. Cole](https://www.phy.cam.ac.uk/directory/colej).

## Citation

If you use **rdfpy** in your work, please cite:

```
@software{rdfpy,
  author       = {Batuhan Yildirim and
                  Hamish Galloway Brown},
  title        = {by256/rdfpy: rdfpy-v0.1.7},
  month        = mar,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v0.1.7},
  doi          = {10.5281/zenodo.4603002},
  url          = {https://doi.org/10.5281/zenodo.4603002}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4298486.svg)](https://doi.org/10.5281/zenodo.4298486)

## Funding

This project was financially supported by the [Science and Technology Facilities Council (STFC)](https://stfc.ukri.org/) and the [Royal Academy of Engineering](https://www.raeng.org.uk/) (RCSRF1819\7\10).

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)




