"""
rdfpy.tests.test_regression
Pins g(r) on the bundled datasets so that future refactors cannot silently
change published numbers.

The reference arrays in ``tests/data/reference/`` were regenerated for the 2.0
release: the values produced by 1.x were wrong for these datasets (both have a
bounding box whose corner is away from the origin, which the 1.x edge
correction mishandled), so they were not preserved.

Regenerate with::

    python -m tests.regenerate_reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
import os
import unittest

import numpy as np

from rdfpy import rdf
from tests.fixtures import brute_force_rdf

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
REFERENCE_DIR = os.path.join(DATA_DIR, 'reference')

# (dataset filename, dr)
DATASETS = [('2dcoords.npy', 5 / 6), ('3dcoords.npy', 0.1)]


def reference_paths(name):
    stem = os.path.splitext(name)[0]
    return (os.path.join(REFERENCE_DIR, stem + '_g_r.npy'),
            os.path.join(REFERENCE_DIR, stem + '_radii.npy'))


class TestBundledDataRegression(unittest.TestCase):

    def test_matches_stored_reference(self):
        for name, dr in DATASETS:
            with self.subTest(dataset=name):
                coords = np.load(os.path.join(DATA_DIR, name))
                g_path, r_path = reference_paths(name)
                self.assertTrue(os.path.exists(g_path),
                                'missing reference for {}; run '
                                'python -m tests.regenerate_reference'
                                .format(name))
                g_r, radii = rdf(coords, dr=dr)
                np.testing.assert_allclose(g_r, np.load(g_path), rtol=1e-10,
                                           equal_nan=True)
                np.testing.assert_allclose(radii, np.load(r_path), rtol=1e-12)

    def test_bundled_data_is_translation_invariant(self):
        """Both bundled datasets sit away from the origin -- 2dcoords starts at
        ~33.9 -- which is exactly the case the 1.x edge correction got wrong."""
        for name, dr in DATASETS:
            with self.subTest(dataset=name):
                coords = np.load(os.path.join(DATA_DIR, name))
                self.assertTrue((coords.min(axis=0) != 0).any(),
                                'dataset no longer exercises the offset case')
                g_ref, _ = rdf(coords, dr=dr)
                g_shift, _ = rdf(coords - coords.min(axis=0), dr=dr)
                np.testing.assert_allclose(g_shift, g_ref, rtol=1e-12,
                                           equal_nan=True)

    def test_2d_bundled_data_matches_brute_force(self):
        """Cross-check the smaller bundled dataset against the independent
        O(N^2) implementation."""
        coords = np.load(os.path.join(DATA_DIR, '2dcoords.npy'))
        g_fast, _ = rdf(coords, dr=5 / 6)
        g_slow, _ = brute_force_rdf(coords, dr=5 / 6)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10, equal_nan=True)


if __name__ == '__main__':
    unittest.main()
