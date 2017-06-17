""" WSFormFactor.py

Woods-Saxon form factor.

Data defined for:
Na-23, I-127, Xe-129, Xe-131, Xe-132, Xe-134, W-184, W-186
arXiv:0608035 (hep-ph)

"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from .FormFactor import FormFactor
from .. import units

import numpy as np

WSdata = {
        (11,23):(2.9393*units.fm,2.994*units.fm,0.523*units.fm),
        (53,127):(5.5931*units.fm,4.749*units.fm,0.523*units.fm),
        (54,129):(5.6315*units.fm,4.776*units.fm,0.523*units.fm),
        (54,131):(5.6384*units.fm,4.781*units.fm,0.523*units.fm),
        (54,132):(5.6460*units.fm,4.787*units.fm,0.523*units.fm),
        (54,134):(5.6539*units.fm,4.792*units.fm,0.523*units.fm),
        (74,184):(6.51*units.fm,5.42*units.fm,0.535*units.fm,
                  0.07*units.fm,0.07*units.fm,0.36*units.fm),#errors
        (74,186):(6.58*units.fm,5.40*units.fm,0.480*units.fm,
                  0.03*units.fm,0.04*units.fm,0.23*units.fm) #errors
       }

class WSFormFactor:
    """ The Woods-Saxon form factor. 

    This is a form factor associated with the charge density
    p(r) = p0 [1 + exp( (r - c)/a )].
 
    There is no closed form solution, so we calculate it
    numerically up to a given tolerance using Romberg
    integration. The Helm form factor ought to be nearly
    identical with the correct choice of parameters.

    """
    def __init__(self, data=WSData[(11,23])):
        """ Initialize with argon as the default. 
 
            Args:
                data: The nuclear data needed for this calculation
                      Default: Sodium-23
                      Given in units of length
        """

        self.c = data[0]
        self.a = data[1]

        self.F0 = 1
        self._tol = 1e-4
    def initialize(self):
        """ Get the normalization constant. """
        self.F0 = 1
        f0 = self.ff2(0)
        self.F0 =  f0

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
                WSData: The data in the correct format and in units
                        of inverse energy.
        """
        if 'WSData' in pars.keys():
            self._data = pars['WSData']
        if 'WSTol' in pars.keys():
            self._tol = pars['WSTol']
        if 'WSa' in pars.keys():
            self.a = pars['WSa']
        if 'WSc' in pars.keys():
            self.c = pars['WSc']

    def _integrand(self,x,q):
        """ The integrand used for integration.

            We transform variables from r to x = arctan(r/a), so our limits
            are 0 to pi/2 and the resulting integrand is finite.
        
            Args:
                x: The unitless parameter arctan(r/a)
                q: The momentum transfer in units of length
        """
        #return 1
        if (x<=0 or x>= 0.5 * np.pi):
            return 0
        r = self.a * np.tan(x)
        #print(r, (r-self.c)/self.a)
        if (q <= 0):
            return (4 * np.pi / self.a * r *r * (1+ (r*r)/(self.a*self.a)) 
                     / (1+np.exp( (r-self.c)/self.a)) ) 


        return (4 * np.pi/q / self.a * r * (1+r*r/self.a/self.a) * np.sin(q*r) 
                 / (1+np.exp( (r-self.c)/self.a)) )

    def ff2(self,Q2):
        """ Calculate the form factor as a function of Q^2. 

            Args:
                Q2: The squared momentum transfer (a positive value)
        """

        q_h = np.sqrt(Q2)/units.hbarc
    
        # Romberg's method of integration
        # We use a tangent change of variables to
        # transform from 0-->infty to 0-->pi/2

        stop_loop = False 
        n = 1
        R = [np.array([0])]
        while(1):
            R.append(np.zeros(n+1) )
            hn = 0.5**(n+1) * np.pi
           
            arr = [self._integrand((2*k-1)*hn,q_h) for k in range(1,2**(n-1)+1)]
            #print(arr)
            
            R[n][0] = 0.5*R[n-1][0] + hn * np.sum(np.array(arr) )
                     
            for m in range(1,n+1):
                R[n][m] = R[n][m-1] + 1./(4**m-1) *(R[n][m-1]-R[n-1][m-1])
        
            if (n > 4
                and np.abs(1-R[n][n]/R[n-1][n-1]) < self._tol
                and np.abs(1-R[n][n]/R[n][n-2]) < self._tol
                and np.abs(1-R[n][n]/R[n][n-3]) < self._tol):
                break
            
            n = n+1
        
        return R[n][n]**2 / self.F0
