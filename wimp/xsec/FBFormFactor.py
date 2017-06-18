""" FBFormFactor.py

Fourier-Bessel form factor. This form factor comes from
a charge density defined by a sum of spherical Bessel
functions j0 with different radial scales.

Data already defined for:
C-12, O-16, Si-28, Si-30, S-32, Ar-40, Ca-40, Ge-70, 
Ge-72, Ge-74, Ge-76, Pb-208
arXiv:0608035
"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from .. import units
from .FormFactor import FormFactor

import numpy as np

_u = units.fm / units.hbarc

FBData = {
    (6,12):(2.464*_u,np.array((0.15721e-1,0.38897e-1,0.37085e-1,0.14795e-1,-0.44831e-2,
            -0.10057e-1,-0.68696e-3,-0.28813e-3,-0.77229e-3,0.66908e-4,
            0.10636e-3,-0.36864e-4,-0.50135e-5,0.94550e-5,-0.47687e-5))*_u,8.0*_u),
    (8,16):(2.737*_u,np.array((0.20238e-1,0.44793e-1,0.33533e-1,0.35030e-2,-0.12293e-1,
            -0.10329e-1,-0.34036e-2,-0.41627e-3,-0.94435e-3,-0.25771e-3,
            0.23759e-3,-0.10603e-3,0.41480e-4))*_u,8.0*_u),
    (14,28):(3.085*_u,np.array((0.33495e-1,0.59533e-1,0.20979e-1,-0.16900e-1,-0.14998e-1,
             -0.93248e-3,0.33266e-2,0.59244e-3,-0.40013e-3,0.12242e-3,
             -0.12994e-4,-0.92784e-5,0.72595e-5,-0.42096e-5))*_u,8.0*_u),
    (14,30):(3.173*_u,np.array((0.28397e-1,0.54163e-1,0.25167e-1,-0.12858e-1,-0.17592e-1,
             -0.46722e-2,0.24804e-2,0.14760e-2,-0.30168e-3,0.483464e-4,
             0.00000,-0.51570e-5,0.30261e-5))*_u,8.5*_u),
    (16,32):(3.248*_u,np.array((0.37251e-1,0.60248e-1,0.14748e-1,-0.18352e-1,-0.10347e-1,
             0.30461e-2,0.35277e-2,-0.39834e-4,-0.97177e-4,0.92279e-4,
             -0.51931e-4,0.22958e-4,-0.86609e-5,0.28879e-5,-0.86632e-6))*_u,8.0*_u),
    (18,40):(3.423*_u,np.array((0.30451e-1,0.55337e-1,0.20203e-1,-0.16765e-1,-0.13578e-1,
             -0.43204e-4,0.91988e-3,-0.41205e-3,0.11971e03,-0.19801e-4,
             -0.43204e-5,0.61205e-5,-0.37803e-5,0.18001e-5,-0.77407e-6))*_u,9.0*_u),
    (20,40):(3.450*_u,np.array((0.44846e-1,0.61326e-1,-0.16818e-2,-0.26217e-1,-0.29725e-2,
             0.85534e-2,0.35322e-2,-0.48259e-3,-0.39346e-3,0.20338e-3,
             0.25461e-4,-0.17794e-4,0.67394e-5,-0.21033e-5))*_u,8.0*_u),
    (32,70):(4.043*_u,np.array((0.38182e-1,0.60306e-1,0.64346e-2,-0.29427e-1,-0.95888e-2,
             0.87849e-2,0.49187e-2,-0.15189e-2,-0.16794e-3,-0.11746e-3,
             0.65768e-4,-0.30691e-4,0.13051e-5,-0.52251e-5))*_u,10.0*_u),
    (32,72):(4.060*_u,np.array((0.38083e-1,0.59342e-1,0.47718e-2,-0.29953e-1,-0.88476e-2,
             0.96205e-2,0.47901e-2,-0.16869e-2,-0.15406e-2,-0.97230e-4,
             -0.47640e-4,-0.15669e-5,0.67076e-5,-0.44500e-5,0.22158e-5))*_u,10.0*_u),
    (32,74):(4.075*_u,np.array((0.37989e-1,0.58298e-1,0.27406e-2,-0.30666e-1,-0.81505e-2,
             0.10231e-1,0.49382e-2,-0.16270e-2,-0.13937e-2,0.15376e-3,
             0.14396e-3,-0.73075e-4,0.31998e-4,-0.12822e-4,0.48406e-5))*_u,10.0*_u),
    (32,76):(4.081*_u,np.array((0.37951e-1,0.57876e-1,0.15303e-2,-0.31822e-1,-0.76875e-2,
             0.11237e-1,0.50780e-2,-0.17293e-2,-0.15523e-2,0.72439e-4,
             0.16560e-3,-0.86631e-4,0.39159e-4,-0.16259e-4,0.63681e-4))*_u,10.0*_u),
    (82,208):(5.499*_u,np.array((0.62732e-1,0.38542e-1,-0.55105e-1,-0.26990e-2,0.31016e-1,
              -0.99486e-2,-0.93012e-2,0.76653e-2,0.20886e-2,-0.17840e-2,
              0.74876e-4,0.32278e-3,-0.11353e-3))*_u,11.0*_u)
         }

def FBFormFactor(FormFactor):
    """ The Fourier-Bessel series form factor.

        This comes from a charge distribution defined by a sum of
        p(r,k) = a(k)j0(k*pi*r / R) over an integer k, with p(r) = 0 
        for r > R.

        Here, j0 is a spherical Bessel function, k goes from 1 to N,
        r is the radius, R is a characteristic size radius, and a(k) 
        are numerical parameters.

        The data should be packaged in such a way that
        data[1] = numpy array of a(k) 
        data[2] = R

        All data should be given in units of inverse energy.

    """
    def __init__(self,data=FBData[(18,40)]):
        """ Initialize with argon as the default. 
 
            Args:
                data: The nuclear data needed for this calculation
                      Default: Argon
                      Given in units of inverse energy (length/ (hbar*c) )

        """
        self._data = data
        self.random = np.random     

    def set_data(self,data):
        """ Set the data.

            Args:
                data: The nuclear data needed for this calculation
                      Given in units of inverse energy (length/ (hbar*c) )
            
        """
        self._data = data

    def set_params(self,pars):
        """ Set the parameters from a dictionray

            Args:
                pars: {string}
 
            Parameters:
                FBData: The data in the correct format and in units
                        of inverse energy.
        """
        if 'FBData' in pars.keys():
            self._data = pars['FBData']

    def ff2(self,Q2):
        """ Calculate the form factor as a function of Q^2. 

            Args:
                Q2: The squared momentum transfer (a positive value)
        """
        if Q2 <=0: 
            return 1
        q = np.sqrt(Q2)

        A = np.sin(q*self._data[2]) / (q*self._data[2])
     
      
        num = np.sum( np.array([ (-1)^(i%2) * self._data[1][i-1] /
                        (i*i*np.pi*np.pi - Q2*self._data[2]**2) 
                        for i in range(1,self._data[1].size+1) ]) )

        denom = np.sum( np.array([ (-1)^(i%2) * self._data[1][i-1] / 
                        (i*i*np.pi*np.pi)
                          for i in range(1,self._data[1].size+1) ]) )

        return A*A *num*num / (denom*denom)

