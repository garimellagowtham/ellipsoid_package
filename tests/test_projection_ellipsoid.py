#!/usr/bin/env python2

import numpy as np
import numpy.testing as np_testing
import unittest
from ellipsoid_package.ellipsoid_projection import *

class ProjectionTest(unittest.TestCase):
    def test_2d_projection(self):
        # Center, principal_axes, radii
        n = 5
        e = Ellipsoid(np.ones(n), np.eye(n), np.arange(n)+1)

        center, Sigma, angle = getXYProjection(n, e)
        np_testing.assert_allclose(center, np.array([1, 1]))
        np_testing.assert_allclose(Sigma, np.array([2, 1]))
        if angle > 180.0:
            angle = angle - 360.0
        elif angle < -180.0:
            angle = angle + 360.0
        np_testing.assert_allclose(np.abs(angle), 90.0)
