#!/usr/bin/env python2
import unittest
import numpy.testing as np_testing
import numpy as np
from ellipsoid_package.ellipsoid_helper import *


class TestEllipsoidProjection(unittest.TestCase):
    def setUp(self):
        np.set_printoptions(precision=3, suppress=True)
        self.n = 3
        
    def testProjectToEllipsoid(self):
        ellipsoid = Ellipsoid(np.zeros(self.n), np.eye(self.n), np.array([4,3,2]))
        x = np.array([1, 0, 0])
        project_x = projectToEllipsoid(x, ellipsoid)
        np_testing.assert_allclose(project_x, np.array([4, 0, 0]))
        x = np.array([0, 0.5, 0])
        project_x = projectToEllipsoid(x, ellipsoid)
        np_testing.assert_allclose(project_x, np.array([0, 1.5, 0]))

    def testProjectToRotatedEllipsoid(self):
        principal_axes = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        ellipsoid = Ellipsoid(np.ones(self.n), principal_axes, np.array([4,3,2]))
        x = np.array([1, 0.4, 0.3])
        project_x = projectToEllipsoid(x, ellipsoid)
        rotated_projected_x = np.dot(principal_axes, project_x  - ellipsoid.center)
        ellipsoid_constraint = np.sum(np.square(rotated_projected_x/ellipsoid.radii))
        self.assertAlmostEqual(ellipsoid_constraint, 1.0, places=3)

    def testProjectSmallVecToEllipsoid(self):
        principal_axes = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        ellipsoid = Ellipsoid(np.ones(self.n), principal_axes, np.array([4,3,2]))
        x = 1e-7*np.ones(3)
        project_x = projectToEllipsoid(x, ellipsoid)
        np_testing.assert_allclose(project_x-ellipsoid.center, np.zeros(3), atol=1e-5)
