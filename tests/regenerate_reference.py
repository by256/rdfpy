"""
rdfpy.tests.regenerate_reference
Regenerates the golden g(r) arrays used by tests/test_regression.py.

Run this only when a change to `rdf` is *intended* to alter published numbers,
and record why in the changelog::

    python -m tests.regenerate_reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
import os

import numpy as np

from rdfpy import rdf
from tests.fixtures import brute_force_rdf
from tests.test_regression import DATA_DIR, DATASETS, REFERENCE_DIR, reference_paths


def main():
    os.makedirs(REFERENCE_DIR, exist_ok=True)
    for name, dr in DATASETS:
        coords = np.load(os.path.join(DATA_DIR, name))
        g_r, radii = rdf(coords, dr=dr)

        # Refuse to enshrine a regression. Without this, a maintainer who
        # breaks `rdf` and then runs the regenerator turns the golden files
        # into a record of the bug and test_matches_stored_reference goes
        # green. Cross-check against the independent O(N^2) reference first.
        g_ref, r_ref = brute_force_rdf(coords, dr=dr)
        try:
            np.testing.assert_allclose(g_r, g_ref, rtol=1e-9, equal_nan=True)
            np.testing.assert_allclose(radii, r_ref, rtol=1e-12)
        except AssertionError as exc:
            raise SystemExit(
                'refusing to write references for {}: rdf disagrees with the '
                'brute-force reference, so the current output is not '
                'trustworthy.\n\n{}'.format(name, exc))

        g_path, r_path = reference_paths(name)
        np.save(g_path, g_r)
        np.save(r_path, radii)
        print('{}: {} bins, g(r) in [{:.4f}, {:.4f}]  (cross-checked against '
              'brute force)'.format(name, len(g_r), np.nanmin(g_r),
                                    np.nanmax(g_r)))


if __name__ == '__main__':
    main()
