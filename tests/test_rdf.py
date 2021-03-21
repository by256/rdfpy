import os
import unittest
import numpy as np
from rdfpy import rdf


class TestRDF2D(unittest.TestCase):

    def test_rdf_2d(self):
        data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/2dcoords.npy')
        coords2d = np.load(data_path)
        # test parallel
        g_r_parallel, radii = rdf(coords2d, dr=5/6, parallel=True)
        # test non-parallel
        g_r, radii = rdf(coords2d, dr=5/6, parallel=False)
        assert np.all(g_r_parallel == g_r)
        

class TestRDF3D(unittest.TestCase):

    def test_rdf_3d(self):
        data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/3dcoords.npy')
        coords3d = np.load(data_path)
        # test parallel
        g_r_parallel, radii = rdf(coords3d, dr=0.1, parallel=True)
        # test non-parallel
        g_r, radii = rdf(coords3d, dr=0.1, parallel=False)
        assert np.all(g_r_parallel == g_r)
