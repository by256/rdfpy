"""
rdfpy.tests.test_physics
Checks g(r) against values that are known analytically, so that a wrong formula
fails a test rather than merely disagreeing with another copy of itself.
~~~~~~~~~~~~~~~~~~~~~~~~
"""
import unittest

import numpy as np
from scipy.spatial import cKDTree

from rdfpy import rdf
from tests.fixtures import (brute_force_rdf, coordination_number, fcc,
                            hexagonal_monolayer, ideal_gas, jittered_lattice,
                            number_density, simple_cubic, square_lattice_2d)


class TestIdealGas(unittest.TestCase):
    """A uniformly random system has g(r) == 1 at every radius, exactly."""

    def _mean_over_seeds(self, d, box, n_particles, dr, n_seeds=6):
        curves = []
        for seed in range(n_seeds):
            coords = ideal_gas(n_particles=n_particles, box=box, d=d, seed=seed)
            g_r, radii = rdf(coords, dr=dr, rho=n_particles / box**d)
            curves.append(g_r)
        return np.array(curves), radii

    def test_ideal_gas_3d(self):
        curves, radii = self._mean_over_seeds(3, 20.0, 8000, 0.5)
        per_bin = curves.mean(axis=0)
        self.assertAlmostEqual(float(curves.mean()), 1.0, delta=0.01)
        np.testing.assert_allclose(per_bin, 1.0, atol=0.05)

    def test_ideal_gas_2d(self):
        curves, radii = self._mean_over_seeds(2, 100.0, 20000, 1.0)
        per_bin = curves.mean(axis=0)
        self.assertAlmostEqual(float(curves.mean()), 1.0, delta=0.01)
        np.testing.assert_allclose(per_bin, 1.0, atol=0.05)

    def test_ideal_gas_with_inferred_rho(self):
        """Not passing rho must still land on 1.0 -- the bounding-box estimate
        is only slightly biased."""
        coords = ideal_gas(n_particles=8000, box=20.0, d=3, seed=11)
        g_r, _ = rdf(coords, dr=0.5)
        self.assertAlmostEqual(float(np.nanmean(g_r)), 1.0, delta=0.02)


class TestInvariances(unittest.TestCase):
    """g(r) depends only on the relative geometry of the particles."""

    def test_translation_invariance(self):
        """Regression test for the edge-correction frame bug: the mask used to
        be evaluated against the pre-translation bounding box, so an offset
        system returned all-NaN."""
        coords = ideal_gas(n_particles=8000, box=20.0, d=3, seed=3)
        g_ref, r_ref = rdf(coords, dr=0.5)

        for offset in (-50.0, 137.0, 1e4):
            with self.subTest(offset=offset):
                g_shift, r_shift = rdf(coords + offset, dr=0.5)
                self.assertFalse(np.isnan(g_shift).any(),
                                 'offset system produced NaN')
                np.testing.assert_allclose(g_shift, g_ref, rtol=1e-12)
                np.testing.assert_allclose(r_shift, r_ref, rtol=1e-12)

    def test_translation_invariance_2d(self):
        coords = ideal_gas(n_particles=20000, box=100.0, d=2, seed=4)
        g_ref, _ = rdf(coords, dr=1.0)
        g_shift, _ = rdf(coords + 500.0, dr=1.0)
        self.assertFalse(np.isnan(g_shift).any())
        np.testing.assert_allclose(g_shift, g_ref, rtol=1e-12)

    def test_anisotropic_translation(self):
        """A different offset per axis must also leave g(r) unchanged."""
        coords = ideal_gas(n_particles=6000, box=20.0, d=3, seed=5)
        g_ref, _ = rdf(coords, dr=0.5)
        g_shift, _ = rdf(coords + np.array([13.0, -71.0, 250.0]), dr=0.5)
        np.testing.assert_allclose(g_shift, g_ref, rtol=1e-12)

    def test_scale_invariance(self):
        """Scaling coordinates and dr by the same factor leaves g(r) unchanged
        and scales the radii."""
        coords = ideal_gas(n_particles=6000, box=20.0, d=3, seed=6)
        g_ref, r_ref = rdf(coords, dr=0.5)
        for factor in (0.25, 7.0):
            with self.subTest(factor=factor):
                g_scaled, r_scaled = rdf(coords * factor, dr=0.5 * factor)
                np.testing.assert_allclose(g_scaled, g_ref, rtol=1e-9)
                np.testing.assert_allclose(r_scaled, r_ref * factor, rtol=1e-9)

    def test_permutation_invariance(self):
        """Reordering the particles must not change the result."""
        coords = ideal_gas(n_particles=4000, box=20.0, d=3, seed=7)
        g_ref, _ = rdf(coords, dr=0.5)
        shuffled = coords[np.random.default_rng(0).permutation(len(coords))]
        g_shuf, _ = rdf(shuffled, dr=0.5)
        np.testing.assert_allclose(g_shuf, g_ref, rtol=1e-12)


class TestLatticeCoordination(unittest.TestCase):
    """Coordination numbers of a lattice are exact integers.

    Regression test for the bin-tiling bug: shells at exact multiples of dr fell
    into the gap between consecutive bins and were silently dropped, so the
    first shell of a simple-cubic lattice integrated to 0 instead of 6.
    """

    def test_simple_cubic_coordination(self):
        dr = 0.05
        coords = simple_cubic(n=20, a=1.0)
        # rcutoff trimmed to just past the outermost shell asserted on; the bins
        # that remain are bit-identical, the rest were never inspected.
        g_r, radii = rdf(coords, dr=dr, rcutoff=0.3)
        rho = number_density(coords)

        shells = [(1.0, 6), (np.sqrt(2), 12), (np.sqrt(3), 8), (2.0, 6),
                  (np.sqrt(5), 24), (np.sqrt(6), 24)]
        for distance, expected in shells:
            with self.subTest(distance=round(distance, 4)):
                n = coordination_number(g_r, radii, dr, rho,
                                        distance - 0.1, distance + 0.1, d=3)
                self.assertAlmostEqual(n, expected, delta=1e-6)

    def test_fcc_first_shell(self):
        dr = 0.02
        coords = fcc(n=10, a=1.0)
        g_r, radii = rdf(coords, dr=dr, rcutoff=0.2)
        rho = number_density(coords)
        n = coordination_number(g_r, radii, dr, rho, 0.6, 0.8, d=3)
        self.assertAlmostEqual(n, 12, delta=1e-6)

    def test_square_lattice_2d_coordination(self):
        dr = 0.05
        coords = square_lattice_2d(n=60, a=1.0)
        g_r, radii = rdf(coords, dr=dr, rcutoff=0.1)
        rho = number_density(coords)

        for distance, expected in [(1.0, 4), (np.sqrt(2), 4), (2.0, 4)]:
            with self.subTest(distance=round(distance, 4)):
                n = coordination_number(g_r, radii, dr, rho,
                                        distance - 0.1, distance + 0.1, d=2)
                self.assertAlmostEqual(n, expected, delta=1e-6)

    def test_hexagonal_monolayer_coordination(self):
        """A 2D lattice whose primitive vectors are not orthogonal and whose
        shells are not integer multiples of the spacing."""
        dr = 0.02
        coords = hexagonal_monolayer(a=1.0)
        g_r, radii = rdf(coords, dr=dr, rho=2 / np.sqrt(3), rcutoff=0.12)

        shells = [(1.0, 6), (np.sqrt(3), 6), (2.0, 6), (np.sqrt(7), 12),
                  (3.0, 6)]
        for distance, expected in shells:
            with self.subTest(distance=round(distance, 4)):
                n = coordination_number(g_r, radii, dr, 2 / np.sqrt(3),
                                        distance - 0.15, distance + 0.15, d=2)
                self.assertAlmostEqual(n, expected, delta=1e-6)

    def test_hexagonal_monolayer_exact_density(self):
        """The known areal density 2/(sqrt(3)*a^2) must make g(r) settle on 1,
        which the bounding-box estimate does not do for this lattice."""
        coords = hexagonal_monolayer(a=1.0)
        g_r, _ = rdf(coords, dr=0.25, rho=2 / np.sqrt(3))
        self.assertAlmostEqual(float(np.nanmean(g_r[len(g_r) // 2:])), 1.0,
                               delta=0.03)

    def test_lattice_peaks_are_present(self):
        """The shells at integer distances must not be missing outright."""
        dr = 0.05
        g_r, radii = rdf(simple_cubic(n=20, a=1.0), dr=dr, rcutoff=0.35)
        for distance in (1.0, 2.0, 3.0):
            with self.subTest(distance=distance):
                near = np.abs(radii - distance) < dr
                self.assertTrue(np.nansum(g_r[near]) > 0.0,
                                'shell at r={} is absent'.format(distance))

    def test_peak_positions_within_half_bin(self):
        """`radii` are bin centres, so a peak is at most dr/2 from the true
        neighbour distance."""
        dr = 0.02
        g_r, radii = rdf(simple_cubic(n=16, a=1.0), dr=dr, rcutoff=0.3)
        for distance in (1.0, np.sqrt(2), np.sqrt(3), 2.0):
            with self.subTest(distance=round(distance, 4)):
                window = np.abs(radii - distance) < 0.1
                peak = radii[window][np.nanargmax(g_r[window])]
                self.assertLessEqual(abs(peak - distance), dr / 2 + 1e-12)


class TestExclusionAndNormalisation(unittest.TestCase):

    def test_hard_sphere_exclusion(self):
        """A jittered lattice has a guaranteed minimum separation; g(r) must be
        exactly zero below it."""
        a, jitter = 1.0, 0.15
        coords = jittered_lattice(n=18, a=a, jitter=jitter, seed=1)
        # The analytic bound a - 2*jitter*sqrt(3) is loose (0.48 here); use the
        # configuration's actual minimum separation so the asserted region
        # covers the whole exclusion zone rather than half of it.
        min_separation = float(cKDTree(coords).query(coords, k=2)[0][:, 1].min())
        self.assertGreater(min_separation, a - 2 * jitter * np.sqrt(3))

        dr = 0.02
        g_r, radii = rdf(coords, dr=dr, rcutoff=0.1)
        below = radii + dr / 2 <= min_separation
        self.assertTrue(below.any(), 'no bins below the exclusion radius')
        np.testing.assert_array_equal(g_r[below], 0.0)

    def test_explicit_rho_scales_output(self):
        """g(r) is inversely proportional to the density it is normalised by."""
        coords = ideal_gas(n_particles=6000, box=20.0, d=3, seed=8)
        inferred = number_density(coords)
        g_default, _ = rdf(coords, dr=0.5)
        g_double, _ = rdf(coords, dr=0.5, rho=2 * inferred)
        np.testing.assert_allclose(g_double, g_default / 2, rtol=1e-12)

    def test_rcutoff_limits_max_radius(self):
        coords = ideal_gas(n_particles=6000, box=20.0, d=3, seed=9)
        dims = coords.max(axis=0) - coords.min(axis=0)
        for rcutoff in (0.5, 0.9, 1.0):
            with self.subTest(rcutoff=rcutoff):
                dr = 0.5
                _, radii = rdf(coords, dr=dr, rcutoff=rcutoff)
                r_max = (dims.min() / 2) * rcutoff
                self.assertLessEqual(radii[-1] + dr / 2, r_max * (1 + 1e-9))
                # and the cutoff is not over-trimmed: one more shell would not
                # have fitted
                self.assertGreater(radii[-1] + 3 * dr / 2, r_max)


class TestAnisotropicBox(unittest.TestCase):
    """r_max and the edge correction key off min(dims). Every other constructed
    fixture is cubic or near-square, so without this the anisotropic path is
    covered only by the bundled datasets."""

    def test_slab_ideal_gas(self):
        rng = np.random.default_rng(21)
        dims = np.array([5.0, 5.0, 50.0])
        coords = rng.uniform(0, 1, size=(6000, 3)) * dims
        g_r, radii = rdf(coords, dr=0.25, rho=6000 / np.prod(dims))
        self.assertLessEqual(radii[-1] + 0.125, (dims.min() / 2) * 0.9 * (1 + 1e-9))
        self.assertAlmostEqual(float(np.nanmean(g_r)), 1.0, delta=0.05)

    def test_slab_matches_brute_force(self):
        rng = np.random.default_rng(22)
        dims = np.array([4.0, 6.0, 30.0])
        coords = rng.uniform(0, 1, size=(1500, 3)) * dims
        g_fast, r_fast = rdf(coords, dr=0.25)
        g_slow, r_slow = brute_force_rdf(coords, dr=0.25)
        np.testing.assert_allclose(r_fast, r_slow, rtol=1e-12)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10, equal_nan=True)

    def test_anisotropic_2d(self):
        rng = np.random.default_rng(23)
        dims = np.array([8.0, 80.0])
        coords = rng.uniform(0, 1, size=(8000, 2)) * dims
        g_r, _ = rdf(coords, dr=0.4, rho=8000 / np.prod(dims))
        self.assertAlmostEqual(float(np.nanmean(g_r)), 1.0, delta=0.05)


class TestAgainstBruteForce(unittest.TestCase):
    """The strongest guard: an independent O(N^2) implementation that shares no
    code with rdfpy."""

    def test_matches_brute_force_3d(self):
        coords = ideal_gas(n_particles=1500, box=10.0, d=3, seed=12)
        g_fast, r_fast = rdf(coords, dr=0.5)
        g_slow, r_slow = brute_force_rdf(coords, dr=0.5)
        np.testing.assert_allclose(r_fast, r_slow, rtol=1e-12)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10,
                                   equal_nan=True)

    def test_matches_brute_force_2d(self):
        coords = ideal_gas(n_particles=1500, box=50.0, d=2, seed=13)
        g_fast, r_fast = rdf(coords, dr=1.0)
        g_slow, r_slow = brute_force_rdf(coords, dr=1.0)
        np.testing.assert_allclose(r_fast, r_slow, rtol=1e-12)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10,
                                   equal_nan=True)

    def test_matches_brute_force_on_lattice(self):
        """Lattices are where the k-d tree bin boundaries are most fragile."""
        coords = simple_cubic(n=8, a=1.0)
        g_fast, _ = rdf(coords, dr=0.05)
        g_slow, _ = brute_force_rdf(coords, dr=0.05)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10,
                                   equal_nan=True)

    def test_matches_brute_force_when_offset(self):
        coords = ideal_gas(n_particles=1200, box=10.0, d=3, seed=14) + 400.0
        g_fast, _ = rdf(coords, dr=0.5)
        g_slow, _ = brute_force_rdf(coords, dr=0.5)
        np.testing.assert_allclose(g_fast, g_slow, rtol=1e-10,
                                   equal_nan=True)


if __name__ == '__main__':
    unittest.main()
