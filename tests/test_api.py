"""
rdfpy.tests.test_api
Input validation, output contract and deprecation behaviour.
~~~~~~~~~~~~~~~~~~~~
"""
import os
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np

from rdfpy import rdf, rdf2d, rdf3d
from tests.fixtures import ideal_gas, simple_cubic


class TestValidation(unittest.TestCase):

    def setUp(self):
        self.coords = ideal_gas(n_particles=2000, box=10.0, d=3, seed=0)

    def test_bad_inputs_raise_value_error(self):
        rng = np.random.default_rng(0)
        cases = {
            '1-D array': (np.array([1.0, 2.0, 3.0]), {}),
            '3-D array': (rng.random((10, 10, 3)), {}),
            '4-D coords': (rng.random((50, 4)), {}),
            '1-D coords': (rng.random((50, 1)), {}),
            'single particle': (np.array([[1.0, 2.0, 3.0]]), {}),
            'empty array': (np.zeros((0, 3)), {}),
            'contains NaN': (np.r_[self.coords, [[np.nan] * 3]], {}),
            'contains inf': (np.r_[self.coords, [[np.inf] * 3]], {}),
            'coplanar': (np.c_[rng.random((200, 2)) * 10, np.zeros(200)], {}),
            'dr zero': (self.coords, {'dr': 0.0}),
            'dr negative': (self.coords, {'dr': -0.5}),
            'dr exceeds box': (self.coords, {'dr': 50.0}),
            'rcutoff zero': (self.coords, {'rcutoff': 0.0}),
            'rcutoff negative': (self.coords, {'rcutoff': -1.0}),
            'rcutoff above 1': (self.coords, {'rcutoff': 2.0}),
            'rho zero': (self.coords, {'rho': 0.0}),
            'rho negative': (self.coords, {'rho': -1.0}),
            # `not rho > 0` alone lets inf through, which silently produced
            # g(r) == 0 in every bin.
            'rho inf': (self.coords, {'rho': np.inf}),
            'rho nan': (self.coords, {'rho': np.nan}),
            'rho non-scalar': (self.coords, {'rho': np.array([1.0, 2.0])}),
            'complex coords': (self.coords.astype(complex), {}),
            'object dtype': (self.coords.astype(object), {}),
            # would ask for ~9e9 bins, i.e. a 67 GiB allocation
            'dr absurdly small': (self.coords, {'dr': 1e-9}),
        }
        for name, (arr, kwargs) in cases.items():
            with self.subTest(case=name):
                params = {'dr': 0.5}
                params.update(kwargs)
                with self.assertRaises(ValueError):
                    rdf(arr, **params)

    def test_validation_survives_O_flag(self):
        """Validation must not be implemented with `assert`, which python -O
        strips. Under -O, 4-D coordinates used to return a silently wrong array
        computed with the 2-D shell-area formula."""
        script = (
            'import numpy as np\n'
            'from rdfpy import rdf\n'
            'try:\n'
            '    rdf(np.random.rand(50, 4), dr=0.05)\n'
            'except ValueError:\n'
            '    raise SystemExit(0)\n'
            'raise SystemExit(1)\n'
        )
        result = subprocess.run([sys.executable, '-O', '-c', script],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0,
                         'validation was skipped under -O:\n' + result.stderr)

    def test_error_messages_are_informative(self):
        with self.assertRaises(ValueError) as ctx:
            rdf(np.random.rand(50, 4), dr=0.05)
        self.assertIn('2 or 3', str(ctx.exception))

        with self.assertRaises(ValueError) as ctx:
            rdf(self.coords, dr=50.0)
        self.assertIn('dr', str(ctx.exception))


class TestInputHandling(unittest.TestCase):

    def test_accepts_nested_lists(self):
        coords = ideal_gas(n_particles=1000, box=10.0, d=3, seed=1)
        g_array, r_array = rdf(coords, dr=1.0)
        g_list, r_list = rdf(coords.tolist(), dr=1.0)
        np.testing.assert_allclose(g_list, g_array, rtol=1e-12)
        np.testing.assert_allclose(r_list, r_array, rtol=1e-12)

    def test_accepts_integer_dtype(self):
        coords = simple_cubic(n=12, a=1.0).astype(np.int64)
        g_r, radii = rdf(coords, dr=0.5)
        self.assertEqual(g_r.shape, radii.shape)
        self.assertTrue(np.isfinite(g_r).any())

    def test_accepts_float32(self):
        coords = ideal_gas(n_particles=2000, box=10.0, d=3,
                           seed=2).astype(np.float32)
        g_r, _ = rdf(coords, dr=1.0)
        self.assertTrue(np.isfinite(g_r).all())

    def test_input_array_is_not_mutated(self):
        coords = ideal_gas(n_particles=1000, box=10.0, d=3, seed=3) + 25.0
        before = coords.copy()
        rdf(coords, dr=1.0)
        np.testing.assert_array_equal(coords, before)


class TestOutputContract(unittest.TestCase):

    def setUp(self):
        self.coords = ideal_gas(n_particles=4000, box=20.0, d=3, seed=4)
        self.dr = 0.5
        self.g_r, self.radii = rdf(self.coords, dr=self.dr)

    def test_shapes_match(self):
        self.assertEqual(self.g_r.shape, self.radii.shape)
        self.assertEqual(self.g_r.ndim, 1)
        self.assertGreater(len(self.g_r), 0)

    def test_radii_are_bin_centres(self):
        self.assertAlmostEqual(self.radii[0], self.dr / 2, places=12)
        np.testing.assert_allclose(np.diff(self.radii), self.dr, rtol=1e-12)

    def test_radii_within_cutoff(self):
        dims = self.coords.max(axis=0) - self.coords.min(axis=0)
        self.assertLessEqual(self.radii[-1] + self.dr / 2,
                             (dims.min() / 2) * 0.9 * (1 + 1e-9))

    def test_exact_quotient_keeps_final_bin(self):
        """r_max/dr an exact integer must yield that many bins, not one fewer
        because the quotient rounded down."""
        rng = np.random.default_rng(0)
        coords = np.r_[rng.uniform(0.5, 9.5, size=(400, 3)),
                       [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]]]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            _, radii = rdf(coords, dr=1.0, rcutoff=1.0)  # r_max/dr == 5 exactly
        self.assertEqual(len(radii), 5)

        # r_max/dr evaluates to 62.99999999999999 here, so without the rounding
        # slack in _BIN_COUNT_TOL this silently loses its last bin.
        _, radii = rdf(simple_cubic(n=8, a=1.0), dr=0.05)
        self.assertEqual(len(radii), 63)

    def test_g_r_is_non_negative(self):
        finite = self.g_r[np.isfinite(self.g_r)]
        self.assertTrue((finite >= 0).all())

    def test_empty_shell_warns_and_returns_nan(self):
        """With the bounding box anchored to exactly [0, 10]^3 and rcutoff=1.0,
        the outermost shell needs an origin exactly at the box centre. No
        particle sits there, so that bin has no valid origins."""
        rng = np.random.default_rng(5)
        coords = np.r_[rng.uniform(0.5, 9.5, size=(500, 3)),
                       [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]]]
        with self.assertWarns(RuntimeWarning):
            g_r, radii = rdf(coords, dr=1.0, rcutoff=1.0)
        self.assertTrue(np.isnan(g_r[-1]))
        self.assertTrue(np.isfinite(g_r[:-1]).all(),
                        'only the outermost bin should be NaN')

    def test_few_origins_warns(self):
        """A bin computed from one or two origins is not an average over
        origins; g(r) there reflects a couple of specific positions."""
        coords = ideal_gas(n_particles=150, box=10.0, d=3, seed=1)
        with self.assertWarns(RuntimeWarning) as ctx:
            g_r, _ = rdf(coords, dr=1.0)
        self.assertIn('origin particles', str(ctx.warning))
        self.assertTrue(np.isfinite(g_r).any())

    def test_warning_message_is_summarised_not_enumerated(self):
        """A fine dr can leave hundreds of unreliable bins; the warning must not
        list every radius."""
        rng = np.random.default_rng(3)
        shell = rng.normal(size=(4000, 3))
        shell /= np.linalg.norm(shell, axis=1, keepdims=True)
        coords = shell * 10.0
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            rdf(coords, dr=0.02)
        messages = [str(w.message) for w in caught
                    if issubclass(w.category, RuntimeWarning)]
        self.assertTrue(messages, 'expected a RuntimeWarning')
        for message in messages:
            self.assertLess(len(message), 300,
                            'warning enumerates bins instead of summarising')

    def test_empty_shell_warns_in_2d(self):
        rng = np.random.default_rng(9)
        coords = np.r_[rng.uniform(0.5, 9.5, size=(400, 2)),
                       [[0.0, 0.0], [10.0, 10.0]]]
        with self.assertWarns(RuntimeWarning):
            g_r, _ = rdf(coords, dr=1.0, rcutoff=1.0)
        self.assertTrue(np.isnan(g_r[-1]))

    def test_duplicate_particles_do_not_self_count(self):
        """Coincident particles sit at distance 0, which belongs to no shell
        since bins are half-open at the inner edge."""
        base = ideal_gas(n_particles=1000, box=10.0, d=3, seed=12)
        g_single, _ = rdf(base, dr=1.0)
        g_double, _ = rdf(np.repeat(base, 2, axis=0), dr=1.0,
                          rho=1000 / np.prod(base.max(0) - base.min(0)))
        # doubling every particle doubles the neighbour count at every r
        np.testing.assert_allclose(g_double, 2 * g_single, rtol=1e-9)

    def test_no_warning_for_well_posed_input(self):
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            warnings.simplefilter('error', DeprecationWarning)
            rdf(self.coords, dr=self.dr)


class TestDeprecations(unittest.TestCase):

    def setUp(self):
        self.coords = ideal_gas(n_particles=2000, box=10.0, d=3, seed=6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.expected, _ = rdf(self.coords, dr=1.0)

    def test_deprecated_kwargs_warn_but_do_not_change_results(self):
        for kwargs in ({'parallel': True}, {'parallel': False},
                       {'progress': True}, {'eps': 1e-15}):
            with self.subTest(kwargs=kwargs):
                with self.assertWarns(DeprecationWarning):
                    g_r, _ = rdf(self.coords, dr=1.0, **kwargs)
                np.testing.assert_allclose(g_r, self.expected, rtol=1e-12)

    def test_parallel_true_and_false_agree(self):
        """The historical contract: both settings give the same numbers.
        Compared with equal_nan, unlike the original test which used `==` and
        would fail spuriously on NaN bins."""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            g_par, r_par = rdf(self.coords, dr=1.0, parallel=True)
            g_ser, r_ser = rdf(self.coords, dr=1.0, parallel=False)
        self.assertTrue(np.allclose(g_par, g_ser, equal_nan=True))
        np.testing.assert_allclose(r_par, r_ser, rtol=1e-12)

    def test_rdf2d_and_rdf3d_warn(self):
        coords3d = ideal_gas(n_particles=1000, box=10.0, d=3, seed=7)
        coords2d = ideal_gas(n_particles=1000, box=10.0, d=2, seed=8)
        with self.assertWarns(DeprecationWarning):
            rdf3d(coords3d, dr=1.0)
        with self.assertWarns(DeprecationWarning):
            rdf2d(coords2d, dr=1.0)

    def test_aliases_reject_mismatched_dimensionality(self):
        coords3d = ideal_gas(n_particles=1000, box=10.0, d=3, seed=9)
        coords2d = ideal_gas(n_particles=1000, box=10.0, d=2, seed=10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            with self.assertRaises(ValueError):
                rdf2d(coords3d, dr=1.0)
            with self.assertRaises(ValueError):
                rdf3d(coords2d, dr=1.0)

    def test_aliases_accept_positional_args(self):
        """In 1.x `rdf2d` was literally `rdf`, so rdf2d(coords, dr, rho,
        rcutoff) was legal. A **kwargs-only shim would break the one function
        whose purpose is backwards compatibility."""
        coords = ideal_gas(n_particles=1000, box=10.0, d=2, seed=14)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            g_pos, _ = rdf2d(coords, 1.0, None, 0.9)
            g_kw, _ = rdf2d(coords, 1.0, rho=None, rcutoff=0.9)
        np.testing.assert_allclose(g_pos, g_kw, rtol=1e-12)

    def test_aliases_pass_through_kwargs(self):
        coords = ideal_gas(n_particles=1000, box=10.0, d=3, seed=15)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            g_alias, _ = rdf3d(coords, 1.0, rho=2.0, rcutoff=0.5)
        g_direct, _ = rdf(coords, 1.0, rho=2.0, rcutoff=0.5)
        np.testing.assert_allclose(g_alias, g_direct, rtol=1e-12)

    def test_aliases_match_rdf(self):
        coords = ideal_gas(n_particles=1000, box=10.0, d=3, seed=11)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            g_alias, _ = rdf3d(coords, dr=1.0)
        g_direct, _ = rdf(coords, dr=1.0)
        np.testing.assert_allclose(g_alias, g_direct, rtol=1e-12)


class TestNoMultiprocessing(unittest.TestCase):
    """rdf must not spawn processes: under the `spawn` start method used on
    macOS and Windows, a Pool inside a library function fails for any caller
    whose script lacks an `if __name__ == '__main__':` guard."""

    # This MUST be written to a real file and executed as a script.
    # `python -c` sets no `__file__` on __main__, so multiprocessing omits
    # `init_main_from_path` from the child's preparation data and the child
    # never re-imports the parent -- which is precisely the recursion this
    # test exists to detect. Run via -c, the test cannot fail.
    SCRIPT = (
        'import multiprocessing as mp\n'
        'mp.set_start_method("spawn", force=True)\n'
        'import numpy as np\n'
        'from rdfpy import rdf\n'
        'g, r = rdf(np.random.default_rng(0).uniform(0, 10, (2000, 3)),'
        ' dr=1.0)\n'
        'assert np.isfinite(g).any()\n'
        'print("ok")\n'
    )

    def test_rdf_runs_under_spawn_start_method_without_main_guard(self):
        env = dict(os.environ)
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env['PYTHONPATH'] = (repo_root + os.pathsep
                             + env.get('PYTHONPATH', ''))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'unguarded.py')
            with open(path, 'w') as f:
                f.write(self.SCRIPT)
            try:
                result = subprocess.run([sys.executable, path], env=env,
                                        capture_output=True, text=True,
                                        timeout=60)
            except subprocess.TimeoutExpired:
                self.fail('rdf hung under the spawn start method in a script '
                          'with no `if __name__ == "__main__":` guard, which '
                          'is what library-level multiprocessing does there')
        self.assertEqual(result.returncode, 0,
                         'rdf failed under the spawn start method:\n'
                         + result.stderr[-2000:])
        self.assertIn('ok', result.stdout)


if __name__ == '__main__':
    unittest.main()
