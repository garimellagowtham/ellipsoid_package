#!/usr/bin/env python

import numpy as np
from collections import namedtuple


Ellipsoid = namedtuple('Ellipsoid', 'center, principal_axes, radii')

def normalize(v):
    norm=np.linalg.norm(v)
    if norm<1e-6:
        return v
    return v/norm

def getProject1Norm(x):
    """
    Find the length of a point projected on an unit cube
    """
    e = normalize(x)
    ub = np.ones_like(x)
    ub_diff = ub - x
    lb_diff = -ub - x
    i_ub_diff = np.argmin(ub_diff)
    i_lb_diff = np.argmax(lb_diff)
    if ub_diff[i_ub_diff] < -lb_diff[i_lb_diff]:
        x_project_1norm = x + (ub_diff[i_ub_diff]/e[i_ub_diff])*e
    else:
        x_project_1norm = x + (lb_diff[i_lb_diff]/e[i_lb_diff])*e
    return np.linalg.norm(x_project_1norm)


def projectToEllipsoid(x, ellipsoid):
    """
    Project a point in a unit cube to an ellipsoid
    """
    scale = 1.0/getProject1Norm(x)
    rotated_x = np.dot(ellipsoid.principal_axes.T, x*ellipsoid.radii)
    return ellipsoid.center + (rotated_x*scale)
