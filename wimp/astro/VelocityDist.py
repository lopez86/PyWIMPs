import numpy as np
from .. import units

class VelocityDist:
    """ This is the base velocity distribution class.
    This particular version implements a truncated 
    Maxwell-Boltzmann distribution:

    f(v) ~ exp( - |v + vE|^2 / v0^2)

    where f(v) = 0 when |v+vE| > vesc.

    This distribution represents a roughly thermalized distribution
    of non-relativistic WIMPs where WIMPs going too fast have
    long escaped the galaxy. The velocity is given in the lab frame,
    so one needs to calculate the velocity of Earth through the
    WIMP halo to get the proper results. The magnitude is needed
    for the recoil energy spectrum, while the direction is required
    if one wants to calculate WIMP directions.

    """

    def __init__(self):
        self._rand = np.random
        self.v0 = 220 * units.km / units.sec
        self.vesc = 550 * units.km / units.sec
        self.vE = 220 * units.km / units.sec * np.array([0,0,1])
        self.norm = 1.0 / (np.pi * self.v0*self.v0)**1.5
        self.needs_norm = True

    def initialize(self):
        """ Do any initial calculations of parameters."""
        self.normalize()


    @property
    def random(self):
        """ The random number generator. """
        return self._rand

    @random.setter
    def random(self,r):
        """ Set the random number generator.
 
            Args:
                r (Numpy RandomState)
        """
        self._rand = r

    def set_params(self,pars):
        """
        Set the parameters for the velocity model.
        
        Args:
            Dictionary with
            'v0' (float) width (1D rms/sqrt(2))
            'vE' (array(3)) velocity of lab frame through DM halo
            'vesc' galactic escape velocity

        """
        # renormalize:
        renorm = False
        if 'v0' in pars:
            self.v0 = pars['v0']
            renorm = True
        if 'vE' in pars:
            self.vE = pars['vE']
        if 'vesc' in pars:
            self.vesc = pars['vesc'] 
            renorm = True
        if renorm:
            self.normalize()

    def f(self,v):
        """
        Calculate the probability density. Call normalize()
        before calling this to get the proper normalization.
         
        Args:
            v (array(3)): WIMP velocity in lab frame

        Returns: 
            float: probability density
        """
        if self.needs_norm:
            self.needs_norm = False
            self.normalize()

        v2 = v+self.vE
        v2 = v2.dot(v2)
        if v2 >= self.vesc*self.vesc:
            return 0
        return self.norm * np.exp( - v2 / (self.v0*self.v0))

    def f_no_escape(self,v):
        """
        Calculates the WIMP velocity probability density
        function ignoring the escape velocity parameter

        Args:
            v (array(3)): The WIMP velocity in the lab frame

        Returns:
            float: probability density
        """
        if self.needs_norm:
            self.needs_norm = False
            self.normalize()
        v2 = v+self.vE
        v2 = v2.dot(v2)
        return self.norm * np.exp( - v2 / (self.v0*self.v0))


    def normalize(self):
        """
        Calculates the normalization constant.
        """

        # First define the limits to integrate over
        # This is non-optimal, but let's just get a rectangular region
        from math import erf

        r = self.vesc / self.v0
        a = np.pi * self.v0**3 * (np.sqrt(np.pi)*erf(r) - 2*r*np.exp(-r*r) )
        self.norm = 1.0 / a
        
        self.needs_norm = False
