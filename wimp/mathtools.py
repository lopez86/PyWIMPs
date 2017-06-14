""" mathtools.py

    This module will hold any mathematical tools that are useful in
    many different files in this package.
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"


import numpy as np

def get_axes(v,tol=1e-4):
""" Given an initial vector, returns an orthogonal basis based on 
    this vector. This has a tolerance parameter, which is used to 
    represent the minimum allowed vector length. If the length is 
    less than this, returns (ex,ey,ez) axes.

    Args:
        v: The vector (3D numpy array)
        tol: The minimum allowed vector length 

    Returns:
        A set of 3 orthogonal unit vectors to be used as a basis.
"""
     e1 = v
     if e1.dot(e1) < tol*tol:
         e1 = np.array([0,0,1])

     e1 = e1 / np.sqrt(e1.dot(e1))
     e2  = np.array([1,0,0])
     e2 = np.cross(e1,e2)
     if e2.dot(e2) < 1e-8:
         e2 = np.array([0,1,0])
     e2 = e2 / np.sqrt(e2.dot(e2))

     e3 = np.cross(e1,e2)
     e3 = e3 / np.sqrt(e3.dot(e3))

     return (e1,e2,e3)
