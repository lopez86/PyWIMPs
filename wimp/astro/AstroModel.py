""" AstroModel.py

    Container for the astrophysics model being used.
"""
__author__ = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

from . import VelocityDist
from .. import units
import numpy as np


class AstroModel:
    """ Container class to hold astrophysics information.

        Attributes:
            velocity: VelocityDist object
    """

    def __init__(self,norm=False):
        """ Initialize the object with some default values.

            Doesn't automatically normalize the velocity
            distribution, although for the standard distribution
            normalization now uses an analytic function.

            Args:
                norm: If true, normalize the velocity distribution.
        """
        self._rand = np.random
        self.velocity = VelocityDist()
        self._wimp_density = 0.3 * units.GeV / (units.cm**3)
        self._vE = 220 * units.km / units.sec * np.array([0,0,1])
        self._v0 = 220 * units.km / units.sec
        self._vesc = 550 * units.km / units.sec 

        self.velocity.vE = self.vE
        self.velocity.v0 = self.v0
        self.velocity.vesc = self.vesc
        
        if norm:
            self.velocity.normalize()

    def initialize(self):
        """ Do any initial calculations of parameters."""
        self.velocity.initialize()


    @property
    def random(self):
        """ The random number generator. """
        return self._rand

    @random.setter
    def set_random(self,r):
        """ Set the random number generator. 
 
            Args: 
                r: (Numpy RandomState)
        """
        self._rand = r
        self.velocity.set_random(r)

    def set_params(self,pars):
        """ Set the parameters using a dictionary.

            Args: 
                pars: {string}
          
            Parameters:
                vE: Earth velocity 3-vector
                v0: Dispersion velocity
                vesc: Galactic escape velocity at Earth
                rhox: The dark matter mass density at Earth
        """
        if 'vE' in pars.keys():
            self.vE = pars['vE']
        if 'v0' in pars.keys():
            self.v0 = pars['v0']
        if 'vesc' in pars.keys():
            self.vesc = pars['vesc']
        if 'rhox' in pars.keys():
            self.wimp_density = pars['rhox']
        self.velocity.set_params(pars)
        
    @property
    def vE(self):
        """ The Earth velocity 3-vector."""
        return self._vE
 
    @vE.setter
    def set_vE(self,ve):
        """ Set the Earth velocity 3-vector.
    
            Args:
                ve: The new vector (numpy.array(3))
        """
        self._vE = ve
        self.velocity.vE = ve

    @property
    def v0(self):
        """ The dispersion velocity. """
        return self._v0

    @v0.setter
    def set_v0(self,v):
        """ Set the dispersion velocity.
  
            Args:
                v: The new velocity
        """
        self._v0 = v
        self.velocity.v0 = v

    @property
    def vesc(self):
        """ The galactic escape velocity. """
        return self._vesc

    @vesc.setter
    def set_vesc(self,v):
        """ Set the galactic escape velocity.
  
            Args:
                v: The new velocity
        """
        self._vesc = v
        self.velocity.vesc = v

    @property
    def wimp_density(self):
        """ The WIMP mass density. """
        return self._wimp_density

    @wimp_density.setter
    def set_wimp_density(self,rho):
        """ Set the WIMP mass density at Earth.
  
            Args:
                v: The new density
        """
        self._wimp_density = rho

    @property
    def rho(self):
        """ The WIMP mass density. """
        return self._wimp_density

