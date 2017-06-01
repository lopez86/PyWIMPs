import numpy as np

def get_axes(v):

     e1 = vE
     if e1.dot(e1) < 1e-8:
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
