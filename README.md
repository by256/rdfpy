# <img alt="rdfpy" src="./logo.png" height="100">

[![Tests](https://github.com/by256/rdfpy/actions/workflows/tests.yml/badge.svg)](https://github.com/by256/rdfpy/actions/workflows/tests.yml)
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

`radii` are **bin centres**: element `k` describes the shell `(k*dr, (k+1)*dr]`, so `radii[0] == dr/2`.

2D coordinates work the same way (the shells become annuli):

```python
coords_2d = np.random.uniform(0.0, 100.0, size=(20000, 2))
g_r, radii = rdf(coords_2d, dr=1.0)
```

**Pass `rho` when you know it.** By default the number density is inferred from the bounding box of
your particles, which stops *at* the outermost ones and so runs high. For disordered systems the bias
is negligible once N is large, but for a finite lattice it is not: a 20×20×20 simple-cubic lattice is
inferred 16.6% too dense, pulling all of `g(r)` down by 14%.

```python
g_r, radii = rdf(coords, dr=0.1, rho=n_particles / box_volume)
```

## Examples

Runnable scripts live in [`examples/`](./examples) (they need matplotlib: `pip install rdfpy[examples]`):

| script | what it shows |
| --- | --- |
| [`monolayer_2d.py`](./examples/monolayer_2d.py) | 2D hexagonal monolayer; generates its own coordinates, and recovers the coordination numbers 6, 6, 6, 12, 6 to within 1% |
| [`crystal.py`](./examples/crystal.py) | 3D crystalline Ti, sharp coordination shells |
| [`water.py`](./examples/water.py) | 3D molecular liquid, decaying structure |

The [Documentation](https://rdfpy.readthedocs.io/) covers choosing `dr`, supplying `rho`, and how
results compare with periodic-boundary tools such as LAMMPS, VMD and freud.

**Note:** In order for **rdfpy** to work correctly, your particles should spatially be in a cuboidal box, where the entire box is filled with particles. Periodic boundary conditions are *not* assumed. Instead, at each radius `r` only those particles whose full shell of radius `r + dr` lies inside the box are used as origins, while every particle remains available as a neighbour.

## How does it work?

- **Fast nearest-neighbor look-up**: a k-d tree is used to count the number of particles as a function of distance from an origin particle.
- **Vectorized, threaded counting**: all origin particles are counted in a single call per shell, and SciPy releases the GIL to spread that work across every available core. No `multiprocessing`, so `rdf` is safe to call from any script on Linux, macOS and Windows without an `if __name__ == "__main__":` guard.

## Upgrading from 1.x

**Version 2.0 changes the numbers `rdf` returns.** Two correctness bugs were fixed and the binning convention was made explicit; if you have published results from 1.x, recompute them.

- **Edge correction is now applied in the right coordinate frame.** In 1.x the interior-particle mask was compared against the pre-translation bounding box. Any system not already anchored at the origin got a wrong `g(r)` and a sufficiently large offset returned all `NaN`.
- **Histogram bins now tile radius space without gaps.** In 1.x each bin covered `(r, r+dr-eps]` while the next started above `r+dr`, so pair distances landing exactly on a bin edge were dropped. This deleted entire coordination shells for crystalline inputs (the first shell of a simple cubic lattice integrated to 0 instead of 6).
- **`radii` are now bin centres starting at `dr/2`**, and the innermost shell `(0, dr]` is now included. In 1.x `radii` held shell left edges starting at `dr`, which shifted every feature down by up to `dr` and never sampled distances below `dr`.
- **Invalid input now raises `ValueError`** instead of `AssertionError` (which `python -O` stripped) or silently returning empty arrays.
- **`parallel`, `progress` and `eps` are deprecated no-ops**; `rdf2d`/`rdf3d` are deprecated in favour of `rdf` and now reject mismatched dimensionality.

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
  note         = {Please cite the version you actually used; results changed in v2.0.0.},
  doi          = {10.5281/zenodo.4625675},
  url          = {https://doi.org/10.5281/zenodo.4625675}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625675.svg)](https://doi.org/10.5281/zenodo.4625675)

## Funding

This project was financially supported by the [Science and Technology Facilities Council (STFC)](https://stfc.ukri.org/) and the [Royal Academy of Engineering](https://www.raeng.org.uk/) (RCSRF1819\7\10).
