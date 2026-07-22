# Changelog

## 2.0.0

**This release changes the values `rdf` returns.** If you have published results produced with 1.x,
recompute them. Two correctness bugs were fixed and the binning convention was made explicit.

### Fixed

- **Edge correction was evaluated in the wrong coordinate frame.** The particles were translated so
  the corner of their bounding box sat at the origin, but the interior-particle mask was still
  compared against the *pre-translation* bounds. Any system not already anchored at the origin got a
  wrong `g(r)`, and a large enough offset returned all `NaN` with only a `RuntimeWarning`. Both
  datasets bundled with the project were affected (up to 69% error in 2D, 93% in 3D). `g(r)` is now
  translation invariant, which is covered by a regression test.
- **Histogram bins did not tile radius space.** Each bin counted the interval `(r, r+dr-eps]` while
  the next began above `r+dr`, leaving the sliver `(r+dr-eps, r+dr]` in no bin at all. Pair distances
  landing exactly on a bin edge were silently discarded: measure-zero for continuum data, but for a
  lattice it deleted entire coordination shells. A simple cubic lattice's first shell integrated to 0
  instead of 6, and 7.3% of the pair distances within the computed radius range were
  discarded outright. Bin edges are now shared between adjacent
  shells and the `eps` nudge is gone.
- **Invalid input is now rejected.** Validation used `assert`, which `python -O` strips: under `-O`,
  4-dimensional coordinates returned a silently wrong array computed with the 2D shell-area formula.
  All validation now raises `ValueError`, and `dr`, `rho`, `rcutoff`, degenerate boxes, non-finite
  coordinates and fewer than two particles are checked too. Several of these previously returned two
  empty arrays with no warning.
- **Divide-by-zero on empty shells.** Shells with no valid origin particles now yield `NaN` and a
  single summarised `RuntimeWarning`, rather than a bare NumPy warning. Shells computed from fewer
  than four origins also warn: averaging over one or two origins is not an average, and the result
  reads as a plausible density while being nothing of the sort.
- **Non-finite `rho` is rejected.** `rho=inf` passed the old `not rho > 0` check and silently
  produced `g(r) == 0` in every bin, while `rho=nan` was rejected (inconsistent, and `inf` is
  exactly what an upstream `N/0` gives you).
- **An absurdly small `dr` now raises** instead of attempting the allocation: `dr=1e-9` in a
  20-unit box asks for 9e9 shells, i.e. a 67 GiB array.
- **Complex and object coordinate dtypes are rejected** rather than silently coerced.

### Changed

- **`radii` are now bin centres**, starting at `dr/2`, and the innermost shell `(0, dr]` is included.
  Previously `radii` held shell left edges starting at `dr`, which shifted every feature down by up to
  `dr` and never sampled distances below `dr`.
- **`multiprocessing` was removed.** Neighbour counting is vectorised over origin particles and
  threaded internally by SciPy. Origins are sorted by their distance to the nearest box face, which
  makes each shell's valid set a prefix of one array (so selecting it is a `searchsorted` rather than
  a fresh `(N, d)` mask per shell) and lets each shell's outer counts serve as the next shell's inner
  counts — one k-d tree query per shell instead of two.

  Measured 1.3–4.7× faster than the old process pool on a 24-thread machine, across
  N = 5 000–100 000 and dr = 0.02–1.0. The gain depends strongly on the number of shells: roughly
  1.3–1.6× at fine `dr` (hundreds of shells) rising to 3–4.7× at coarse `dr`.

  This also fixes `rdf` being unusable under the `spawn` start method (the default on macOS and
  Windows) from any script lacking an `if __name__ == "__main__":` guard.
- `matplotlib` moved from `install_requires` to `extras_require['examples']`; nothing in the package
  imports it.
- `find_packages()` now excludes `tests`, which was previously installed as a top-level package into
  `site-packages`.
- Source distributions now include `tests/` and its data, so the suite can be run from an sdist.
- CI moved from Travis (shut down in 2021, so the project effectively had none) to GitHub Actions,
  covering Linux, macOS and Windows on Python 3.9–3.13, including a `python -O` run.

### Deprecated

- `parallel`, `progress` and `eps` are accepted but ignored, and emit `DeprecationWarning`.
- `rdf2d` and `rdf3d` emit `DeprecationWarning` and now reject coordinates of the wrong
  dimensionality. Previously they were plain aliases of `rdf`, so `rdf2d` silently accepted 3D
  input. They still accept positional arguments, since in 1.x they *were* `rdf`.

### Documentation and examples

- Added a **2D example** (`examples/monolayer_2d.py`), the first anywhere in the project despite 2D
  being a headline feature. It builds a hexagonal monolayer in-script rather than shipping another
  `.npy`, and verifies itself by recovering the coordination numbers 6, 6, 6, 12, 6 to within 1%.
- Moved `rdfpy/examples/` to a top-level `examples/`. It sat inside the package directory but had no
  `__init__.py`, so it shipped in neither the wheel nor the sdist: ~1.4 MB of files in a folder that
  looked installed but was not.
- Regenerated the example figures, which had been produced by an older revision and no longer matched
  what the code output.
- New documentation on choosing `dr` (with measured resolution/noise figures), when to supply `rho`,
  and how results differ from periodic-boundary tools (LAMMPS, VMD, MDAnalysis, freud): the most
  likely source of "why don't my numbers match".
- Corrected the documented formula: it gave the 3D-only thin-shell approximation `4*pi*r^2*dr` while
  the code uses exact shell volumes, omitted the average over origin particles, never mentioned the
  edge-correction scheme, and wrote "between distances r and dr" for "r and r+dr".

### Testing

The 1.x suite consisted of two tests that only asserted the parallel and serial paths agreed. Since
both ran the same code, they passed against every bug above. The suite is now 54 tests pinning `g(r)`
to analytically known values: ideal-gas convergence to 1.0 in 2D and 3D, translation/scale/permutation
invariance, exact coordination numbers for simple cubic (6, 12, 8, 6, 24, 24), FCC (12), 2D square (4, 4, 4) and 2D
hexagonal (6, 6, 6, 12, 6) lattices, hard-sphere exclusion, and agreement with an independent brute-force O(N²)
implementation whose pair counting, box translation and shell-count contract are derived
independently of `rdfpy`.

The suite was validated by mutation testing: 34 deliberate bugs were injected into `rdf`, and every
one that changes behaviour is caught. Two gaps that testing exposed have been closed — the
`spawn`-safety test ran its script via `python -c`, which sets no `__file__` on `__main__` and so
structurally exempts the child from the re-import recursion the test exists to detect (it could not
fail); and the bin-count rounding guard was covered only by coincidence in an unrelated fixture.

## 1.0.0

Initial stable release.
