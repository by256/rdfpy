import warnings

import numpy as np
from scipy.spatial import cKDTree


_UNSET = object()

# Slack used when working out how many whole shells fit inside r_max, so that a
# quotient which should be an integer is not lost to floating-point rounding
# (e.g. r_max/dr evaluating to 62.99999999999999 rather than 63).
_BIN_COUNT_TOL = 1e-9

# Below this many origin particles, averaging over origins is not meaningfully
# an average and g(r) reflects a couple of specific positions. Deliberately set
# low: `rcutoff` already trims the merely-noisy outermost bins, so this is here
# to flag the degenerate case, not to editorialise about ordinary noise.
_MIN_VALID_ORIGINS = 4

# Guard against a mistyped dr silently requesting an unallocatable histogram:
# dr=1e-9 in a 20-unit box asks for 9e9 bins, i.e. a 67 GiB array.
_MAX_BINS = 1_000_000


def _deprecated(name, value, replacement=None, stacklevel=3):
    """Warn that a keyword argument is accepted but ignored."""
    if value is _UNSET:
        return
    msg = '`{}` is deprecated and has no effect; it will be removed in a ' \
          'future release.'.format(name)
    if replacement is not None:
        msg += ' ' + replacement
    warnings.warn(msg, DeprecationWarning, stacklevel=stacklevel)


def _validate(particles, dr, rho, rcutoff):
    """Validate and normalize the inputs to `rdf`.

    Returns the particle coordinates as a float64 array, translated so that the
    corner of their bounding box sits at the origin, along with the box
    dimensions and the number density.
    """
    particles = np.asarray(particles)
    if particles.dtype == object or not np.issubdtype(particles.dtype,
                                                      np.number):
        raise ValueError('particles must have a numeric dtype; got dtype '
                         '{}.'.format(particles.dtype))
    if np.issubdtype(particles.dtype, np.complexfloating):
        raise ValueError('particles must be real; got dtype '
                         '{}.'.format(particles.dtype))
    particles = particles.astype(np.float64)

    if particles.ndim != 2:
        raise ValueError(
            'particles must be an array of shape (N, d), where N is the number '
            'of particles and d is the number of dimensions; got an array with '
            '{} dimension(s).'.format(particles.ndim))

    N, d = particles.shape
    if d not in (2, 3):
        raise ValueError('RDF can only be computed in 2 or 3 dimensions; got '
                         '{}-dimensional coordinates.'.format(d))
    if N < 2:
        raise ValueError('at least 2 particles are required to compute an RDF; '
                         'got {}.'.format(N))
    if not np.all(np.isfinite(particles)):
        raise ValueError('particles contains NaN or infinite coordinates.')

    if not np.isscalar(dr) and np.asarray(dr).ndim != 0:
        raise ValueError('dr must be a scalar; got {!r}.'.format(dr))
    if not dr > 0:
        raise ValueError('dr must be greater than 0; got {}.'.format(dr))
    if not 0 < rcutoff <= 1:
        raise ValueError('rcutoff must be in the interval (0, 1]; got '
                         '{}.'.format(rcutoff))

    mins = particles.min(axis=0)
    dims = particles.max(axis=0) - mins
    if not np.all(dims > 0):
        axes = np.flatnonzero(dims <= 0).tolist()
        raise ValueError(
            'particles span zero extent along axis/axes {}, so the box volume '
            'is zero. Drop the degenerate axis and pass {}-dimensional '
            'coordinates instead.'.format(axes, d - len(axes)))

    # Translate so that the box is [0, dims]. Every subsequent geometric test is
    # made in this frame -- mins/maxs must not be used again.
    particles = particles - mins

    if rho is None:
        rho = N / np.prod(dims)
    else:
        rho_arr = np.asarray(rho, dtype=np.float64)
        if rho_arr.ndim != 0:
            raise ValueError('rho must be a scalar; got shape '
                             '{}.'.format(rho_arr.shape))
        rho = float(rho_arr)
        # `not rho > 0` alone would let inf through, which silently yields
        # g(r) == 0 everywhere.
        if not (rho > 0 and np.isfinite(rho)):
            raise ValueError('rho must be a finite number greater than 0; got '
                             '{}.'.format(rho))

    return particles, dims, float(rho)


def rdf(particles, dr, rho=None, rcutoff=0.9, eps=_UNSET, parallel=_UNSET,
        progress=_UNSET):
    """
    Computes the 2D or 3D radial distribution function g(r) of a set of
    particle coordinates of shape (N, d). Particles must be placed in a 2D or 3D
    cuboidal box of dimensions [width x height (x depth)].

    No periodic boundary conditions are assumed. Instead, an edge correction is
    applied: at radius r only those particles whose entire shell of radius
    r + dr lies inside the box are used as origins, while all particles remain
    available as neighbours.

    Parameters
    ----------
    particles : (N, d) np.array
        Set of particles from which to compute the radial distribution function
        g(r). Must be of shape (N, 2) or (N, 3) for 2D and 3D coordinates
        respectively.
    dr : float
        Delta r. Determines the spacing between successive radii over which g(r)
        is computed.
    rho : float, optional
        Number density. If left as None, the box dimensions are inferred from
        the bounding box of the particles and the number density is calculated
        accordingly. That bounding box stops at the outermost particles, so the
        inferred value runs high -- by roughly (1 + 2d/N) for a disordered
        system, but by a whole lattice spacing per axis for a crystal (16.6% for
        a 20x20x20 simple-cubic lattice). Pass rho explicitly whenever it is
        known.
    rcutoff : float
        Radii cutoff value in the interval (0, 1]. g(r) is computed out to
        rcutoff * min(box dimensions) / 2, i.e. the default of 0.9 stops at 90%
        of half the shortest box dimension. This trims the largest radii, which
        are computed from the fewest origin particles and are correspondingly
        noisy.
    eps : deprecated
        No longer used. Bin edges are now exact, so no epsilon nudge is needed.
    parallel : deprecated
        No longer used. Neighbour counting is vectorised and threaded internally
        by SciPy across all available cores.
    progress : deprecated
        No longer used.

    Returns
    -------
    g_r : (n_radii,) np.array
        Radial distribution function values g(r). Bins with no valid origin
        particles are set to NaN and raise a RuntimeWarning; bins computed from
        very few origins also raise a RuntimeWarning.
    radii : (n_radii,) np.array
        Radii over which g(r) is computed. These are **bin centres**: element k
        describes the spherical shell (k*dr, (k+1)*dr], so radii[0] == dr/2.

    Raises
    ------
    ValueError
        If the coordinates, dr, rho or rcutoff are invalid.

    Examples
    --------
    >>> import numpy as np
    >>> from rdfpy import rdf
    >>> coords = np.random.uniform(0.0, 20.0, size=(2500, 3))
    >>> g_r, radii = rdf(coords, dr=0.5)
    """
    _deprecated('eps', eps)
    _deprecated('parallel', parallel,
                'Neighbour counting is now vectorised and threaded internally.')
    _deprecated('progress', progress)

    particles, dims, rho = _validate(particles, dr, rho, rcutoff)
    N, d = particles.shape

    r_max = (np.min(dims) / 2) * rcutoff
    n_bins = int(np.floor(r_max / dr + _BIN_COUNT_TOL))
    if n_bins < 1:
        raise ValueError(
            'dr ({}) is larger than the largest usable radius ({:.6g} = '
            'rcutoff * min(box dimensions) / 2), so no shells can be computed. '
            'Use a smaller dr.'.format(dr, r_max))
    if n_bins > _MAX_BINS:
        raise ValueError(
            'dr ({}) would require {} shells out to r_max = {:.6g}, which is '
            'above the limit of {}. Use a larger dr.'.format(
                dr, n_bins, r_max, _MAX_BINS))

    # A single edge array guarantees that consecutive shells (edges[k],
    # edges[k+1]] tile radius space exactly: no gaps, no overlaps. Sharing the
    # endpoints is what keeps pair distances that land exactly on a bin
    # boundary -- every coordination shell of a lattice -- from being dropped.
    edges = np.arange(n_bins + 1) * dr
    radii = edges[:-1] + dr / 2

    tree = cKDTree(particles)

    # An origin is valid for a shell when its whole shell fits inside the box,
    # i.e. when its distance to the nearest face is at least the outer radius.
    # Sorting by that margin makes the valid sets nested prefixes of one array,
    # so selecting them is a searchsorted rather than a fresh (N, d) mask per
    # bin, and each bin's outer counts can be reused as the next bin's inner
    # counts -- one tree query per bin instead of two.
    margin = np.minimum(particles, dims - particles).min(axis=1)
    order = np.argsort(-margin, kind='stable')
    sorted_particles = particles[order]
    neg_sorted_margin = -margin[order]          # ascending, for searchsorted

    g_r = np.empty(n_bins)
    empty_bins = []
    sparse_bins = []
    prev_counts = None
    prev_n_valid = -1

    for k in range(n_bins):
        r_inner, r_outer = edges[k], edges[k + 1]

        n_valid = int(np.searchsorted(neg_sorted_margin, -r_outer,
                                      side='right'))
        if n_valid == 0:
            g_r[k] = np.nan
            empty_bins.append(k)
            prev_counts, prev_n_valid = None, 0
            continue
        if n_valid < _MIN_VALID_ORIGINS:
            sparse_bins.append(k)

        valid_particles = sorted_particles[:n_valid]

        if prev_counts is not None and prev_n_valid >= n_valid:
            # The previous bin's outer edge is this bin's inner edge, and its
            # valid set is a superset of this one's, so its counts still apply.
            inner_counts = prev_counts[:n_valid]
        else:
            inner_counts = tree.query_ball_point(
                valid_particles, r_inner, return_length=True, workers=-1)
        outer_counts = tree.query_ball_point(
            valid_particles, r_outer, return_length=True, workers=-1)

        # Each origin counts itself in both terms, so self-counting cancels.
        n = (outer_counts - inner_counts).sum()

        if d == 3:
            shell_vol = (4 / 3) * np.pi * (r_outer**3 - r_inner**3)
        else:
            shell_vol = np.pi * (r_outer**2 - r_inner**2)

        g_r[k] = n / (n_valid * shell_vol * rho)
        prev_counts, prev_n_valid = outer_counts, n_valid

    _warn_about_unreliable_bins(empty_bins, sparse_bins, radii, n_bins)

    return g_r, radii


def _warn_about_unreliable_bins(empty_bins, sparse_bins, radii, n_bins):
    """Summarise bins that have no or few origin particles.

    Summarised rather than enumerated: a fine dr can leave hundreds of such
    bins, and a warning listing every radius is unreadable.
    """
    if empty_bins:
        warnings.warn(
            '{} of {} bins (r = {:.4g} to {:.4g}) have no particles far enough '
            'from the box edges to act as origins; g(r) is NaN there. Reduce '
            'rcutoff or use a larger box.'.format(
                len(empty_bins), n_bins, radii[empty_bins[0]],
                radii[empty_bins[-1]]),
            RuntimeWarning, stacklevel=3)
    if sparse_bins:
        warnings.warn(
            '{} of {} bins (r = {:.4g} to {:.4g}) are computed from fewer than '
            '{} origin particles, so g(r) there reflects a handful of specific '
            'positions rather than a density. Reduce rcutoff or use a larger '
            'box.'.format(
                len(sparse_bins), n_bins, radii[sparse_bins[0]],
                radii[sparse_bins[-1]], _MIN_VALID_ORIGINS),
            RuntimeWarning, stacklevel=3)


def rdf2d(particles, dr, *args, **kwargs):
    """Deprecated alias for :func:`rdf` restricted to 2D coordinates."""
    warnings.warn('rdf2d is deprecated; use rdf instead.', DeprecationWarning,
                  stacklevel=2)
    particles = np.asarray(particles)
    if particles.ndim != 2 or particles.shape[-1] != 2:
        raise ValueError('rdf2d requires coordinates of shape (N, 2); got '
                         '{}. Use rdf instead.'.format(particles.shape))
    return rdf(particles, dr, *args, **kwargs)


def rdf3d(particles, dr, *args, **kwargs):
    """Deprecated alias for :func:`rdf` restricted to 3D coordinates."""
    warnings.warn('rdf3d is deprecated; use rdf instead.', DeprecationWarning,
                  stacklevel=2)
    particles = np.asarray(particles)
    if particles.ndim != 2 or particles.shape[-1] != 3:
        raise ValueError('rdf3d requires coordinates of shape (N, 3); got '
                         '{}. Use rdf instead.'.format(particles.shape))
    return rdf(particles, dr, *args, **kwargs)
