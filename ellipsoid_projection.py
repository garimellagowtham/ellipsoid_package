#!/usr/bin/env python2
from matplotlib.patches import Ellipse
from ellipsoid_tool import EllipsoidTool
from ellipsoid_helper import Ellipsoid
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
    print "Center: ", center
    print "Sigma: ", Sigma
    print "angle: ", angle
    e = Ellipse(center, 2*Sigma[0], 2*Sigma[1], angle, fill=False)
    subplot.add_patch(e)

def plotEllipsoid3D(n, ellipse, subplot):
    center, Sigma, rotation = getXYZProjection(n, ellipse)
    print "Center: ", center
    print "Sigma: ", Sigma
    #print "rotation: ", rotation
    ellipsoid_tool.plotEllipsoid(center, Sigma, rotation, subplot, cageColor='k')
