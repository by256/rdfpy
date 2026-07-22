"""
rdfpy.tests.fixtures
Particle configurations with analytically known radial distribution functions.
~~~~~~~~~~~~~~~~~~~~
"""
import numpy as np


def simple_cubic(n=20, a=1.0):
    """Simple-cubic lattice, n**3 sites, spacing `a`.

    Coordination shells (distance, count): (a, 6), (a*sqrt2, 12), (a*sqrt3, 8),
    (2a, 6), (a*sqrt5, 24), (a*sqrt6, 24).
    """
    grid = np.arange(n) * a
    x, y, z = np.meshgrid(grid, grid, grid, indexing='ij')
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)


def fcc(n=10, a=1.0):
    """Face-centred-cubic lattice, 4*n**3 sites, conventional cell length `a`.

    First coordination shell: 12 neighbours at a/sqrt(2).
    """
    basis = np.array([[0.0, 0.0, 0.0],
                      [0.0, 0.5, 0.5],
                      [0.5, 0.0, 0.5],
                      [0.5, 0.5, 0.0]])
    grid = np.arange(n)
    i, j, k = np.meshgrid(grid, grid, grid, indexing='ij')
    cells = np.stack([i.ravel(), j.ravel(), k.ravel()], axis=1)
    return ((cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)) * a


def square_lattice_2d(n=60, a=1.0):
    """2D square lattice, n**2 sites, spacing `a`.

    Coordination shells: (a, 4), (a*sqrt2, 4), (2a, 4), (a*sqrt5, 8).
    """
    grid = np.arange(n) * a
    x, y = np.meshgrid(grid, grid, indexing='ij')
    return np.stack([x.ravel(), y.ravel()], axis=1)


def hexagonal_monolayer(n_x=60, n_y=70, a=1.0):
    """2D triangular (hexagonal) lattice, rows offset by a/2.

    Unlike the square lattice its primitive vectors are not orthogonal, so it
    exercises the 2D normalisation on a structure whose shells are not at
    integer multiples of the spacing. Coordination shells (distance, count):
    (a, 6), (a*sqrt3, 6), (2a, 6), (a*sqrt7, 12), (3a, 6). Areal density is
    exactly 2/(sqrt(3)*a**2).
    """
    i, j = np.meshgrid(np.arange(n_x), np.arange(n_y), indexing='ij')
    x = (i + 0.5 * (j % 2)) * a
    y = j * (np.sqrt(3) / 2) * a
    return np.stack([x.ravel(), y.ravel()], axis=1)


def ideal_gas(n_particles=8000, box=20.0, d=3, seed=0):
    """Uniformly random ("ideal gas") coordinates. g(r) == 1 at every radius."""
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, box, size=(n_particles, d))


def jittered_lattice(n=18, a=1.0, jitter=0.15, seed=0):
    """Simple-cubic lattice with bounded uniform displacements.

    Every pair separation is at least ``a - 2*jitter*sqrt(3)``, so g(r) must be
    exactly zero below that distance.
    """
    rng = np.random.default_rng(seed)
    sites = simple_cubic(n=n, a=a)
    return sites + rng.uniform(-jitter, jitter, size=sites.shape)


def shell_volumes(radii, dr, d):
    """Volume (3D) or area (2D) of the shell each bin centre describes."""
    inner = radii - dr / 2
    outer = radii + dr / 2
    if d == 3:
        return (4 / 3) * np.pi * (outer**3 - inner**3)
    return np.pi * (outer**2 - inner**2)


def coordination_number(g_r, radii, dr, rho, lo, hi, d=3):
    """Integrate rho * g(r) * V_shell(r) over the radii in [lo, hi).

    For a lattice this returns the number of neighbours in that coordination
    shell, which is an exact integer known from the crystal structure.
    """
    mask = (radii >= lo) & (radii < hi)
    return float(np.nansum((rho * g_r * shell_volumes(radii, dr, d))[mask]))


def number_density(particles):
    """Number density implied by the bounding box, as `rdf` infers it."""
    particles = np.asarray(particles, dtype=float)
    dims = particles.max(axis=0) - particles.min(axis=0)
    return particles.shape[0] / np.prod(dims)


def brute_force_rdf(particles, dr, rho=None, rcutoff=0.9):
    """O(N^2) reference RDF, partially independent of ``rdfpy``.

    Be precise about what this does and does not establish. Independently
    derived here, and therefore genuinely cross-checked:

    * the pair counting -- a dense ``np.linalg.norm`` distance matrix rather
      than k-d tree count differencing (this is where bug C2 lived);
    * the translation into the box frame (where bug C1 lived);
    * the bin-count contract, counted up explicitly rather than by division.

    Deliberately restated rather than re-derived, and therefore NOT
    cross-checked by comparison: the edge-correction mask, the shell volume,
    and the normalisation. An error in one of those *contracts* would change
    both implementations together and this comparison would stay green --
    ``TestIdealGas`` and the stored references in ``test_regression`` are what
    catch that class of error.
    """
    particles = np.asarray(particles, dtype=float)
    mins = particles.min(axis=0)
    dims = particles.max(axis=0) - mins
    particles = particles - mins
    n, d = particles.shape
    if rho is None:
        rho = n / np.prod(dims)

    # Same contract as rdf -- the number of whole shells of width dr that fit
    # inside r_max -- but counted up explicitly rather than by division, so a
    # quotient like 62.99999999999999 does not cost a bin.
    r_max = (dims.min() / 2) * rcutoff
    n_bins = int(np.floor(r_max / dr))
    while (n_bins + 1) * dr <= r_max * (1 + 1e-9):
        n_bins += 1
    edges = np.arange(n_bins + 1) * dr

    distances = np.linalg.norm(particles[:, None, :] - particles[None, :, :],
                               axis=2)

    g_r = np.empty(n_bins)
    for k in range(n_bins):
        r_inner, r_outer = edges[k], edges[k + 1]
        valid = np.all((particles >= r_outer) & (particles <= dims - r_outer),
                       axis=1)
        if not valid.any():
            g_r[k] = np.nan
            continue
        sub = distances[valid]
        count = ((sub > r_inner) & (sub <= r_outer)).sum()
        if d == 3:
            shell = (4 / 3) * np.pi * (r_outer**3 - r_inner**3)
        else:
            shell = np.pi * (r_outer**2 - r_inner**2)
        g_r[k] = count / (valid.sum() * shell * rho)

    return g_r, edges[:-1] + dr / 2
