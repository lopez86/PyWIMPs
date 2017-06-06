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
        self.v0 = 220 * units.km / units.sec
        self.vesc = 550 * units.km / units.sec
        self.vE = 220 * units.km / units.sec * np.array([0,0,1])
        self.norm = 1.0 / (np.pi * self.v0*self.v0)**1.5
        self.max_iter = 1000000 # Gets ~0.01% error for fairly standard assumptions
        self.tol_norm = 0.001

    def set_params(self,pars):
        """
        Set the parameters for the velocity model.
        
        Args:
            Dictionary with
            'v0' (float) width (1D rms/sqrt(2))
            'vE' (array(3)) velocity of lab frame through DM halo
            'vesc' galactic escape velocity
        """
        if 'v0' in pars.keys():
            self.v0 = pars['v0']
        if 'vE' in pars.keys():
            self.vE = pars['vE']
        if 'vesc' in pars.keys():
            self.vesc = pars['vesc']
        self.norm = 1.0 / (np.pi * self.v0*self.v0)**1.5
      

    def f(self,v):
        """
        Calculate the probability density. Call normalize()
        before calling this to get the proper normalization.
         
        Args:
            v (array(3)): WIMP velocity in lab frame

        Returns: 
            float: probability density
        """
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
        v2 = v+self.vE
        v2 = v2.dot(v2)
        return self.norm * np.exp( - v2 / (self.v0*self.v0))


    def normalize(self,calcErr=False):
        """
        Monte Carlo integration of f(v). Numerically calculates the
        integral of f(v) and then renormalizes f(v) so that the 
        integral is set to unity.

        Args:
            calcErr (bool): If True, calculate the approximate error

        """
        # First define the limits to integrate over
        # This is non-optimal, but let's just get a rectangular region
        
        niter = 0
        fave = 0

        while niter < self.max_iter:
            if niter%(self.max_iter/10) == 0:
                print("Normalization: %0.1f%% done"%(100 * (niter / (self.max_iter)) ,))
            vec = np.random.normal(-self.vE,self.v0/np.sqrt(2),3)
            vec2 = vec + self.vE
            vec_prob = 1./(np.pi*self.v0*self.v0)**1.5 * \
                       np.exp( - (vec2.dot(vec2)) / (self.v0*self.v0) )
            fval = self.f(vec) / vec_prob
            

            fave = fave + fval

            niter = niter + 1
        
        fave = fave / niter


       
        if calcErr is True:
            fvar = 0
            niter = 0
            while niter < self.max_iter:
                vec = np.random.normal(-self.vE,self.v0/np.sqrt(2),3)
                vec2 = vec + self.vE
                vec_prob = 1./(np.pi*self.v0*self.v0)**1.5 * \
                       np.exp( - (vec2.dot(vec2)) / (self.v0*self.v0) )
            
                fvar = fvar + (self.f(vec)/vec_prob-fave)**2
                niter = niter + 1
            fvar = fvar / (self.max_iter*(self.max_iter-1))
            ferr = np.sqrt(fvar) / fave
            print("Norm: %f Err: %f"%(fave,ferr))
            print("Fractional normalization error is:")
            print(ferr)
            if ( ferr > self.tol_norm):
                print(("This is above tolerance. Please increase the number of "
                       "iterations"))

        self.norm = self.norm / (fave)
