#!/usr/bin/env python2
from matplotlib.patches import Ellipse
from ellipsoid_tool import EllipsoidTool
import numpy as np


ellipsoid_tool = EllipsoidTool()

def getXYProjection(n, ellipse):
    #Projection matrix
    e1 = np.zeros(n)
    e1[0] = 1
    e2 = np.zeros(n)
    e2[1] = 1
    T = np.vstack((e1, e2))
    L = np.dot(ellipse.principal_axes.T, np.diag(ellipse.radii))
    U,Sigma,V = np.linalg.svd(np.dot(T, L))
    if np.linalg.det(U) < 0:
        U[:,0] = -1*U[:,0]
    angle = np.arctan2(U[1, 0], U[0, 0])
    center = np.dot(T, ellipse.center)
    return (center, Sigma, angle*180.0/np.pi)

def getXYZProjection(n, ellipse):
    #Projection matrix
    e1 = np.zeros(n)
    e1[0] = 1
    e2 = np.zeros(n)
    e2[1] = 1
    e3 = np.zeros(n)
    e3[2] = 1
    T = np.vstack((e1, e2, e3))
    L = np.dot(ellipse.principal_axes.T, np.diag(ellipse.radii))
    U,Sigma,V = np.linalg.svd(np.dot(T, L))
    if np.linalg.det(U) < 0:
        U[:,0] = -1*U[:,0]
    center = np.dot(T, ellipse.center)
    return (center, Sigma, U)

def plotEllipse(n, ellipse, subplot):
    center, Sigma, angle = getXYProjection(n, ellipse)
    e = Ellipse(center, 2*Sigma[0], 2*Sigma[1], angle, fill=False,
                ec='k')
    subplot.add_patch(e)

def plotEllipsoid3D(n, ellipse, subplot):
    center, Sigma, rotation = getXYZProjection(n, ellipse)
    ellipsoid_tool.plotEllipsoid(center, Sigma, rotation, subplot, cageColor='k')

# Obtained from stack overflow:
# https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)
