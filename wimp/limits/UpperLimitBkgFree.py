""" UpperLimitBkgFree

    Calculate frequentist upper limits given an all-signal
    (background-free) model.
"""
from .. import units
from ..astro.VelocityDist import VelocityDist
from ..astro.AstroModel import AstroModel
from ..xsec.CrossSection import CrossSection
from ..xsec.FormFactor import FormFactor
from ..xsec.InteractionModel import InteractionModel
from ..mc.MaxwellWeightedSampler import MaxwellWeightedSampler
from ..mc.sample import Sample

import numpy as np
import scipy.stats

class UpperLimitBkgFree:
    """ Class to calculate simple upper limits where it is 
        assumed that there are no background events.

        For counting experiments.

        Attributes:
            _tol: The fractional tolerance for rate calculations 
                  and cross section limits.
            model: The physics model used to calculate rates. 
                   Must have sample() defined.
            _Emin: Minimum energy
            _Emax: Maximum energy 
 
    """
    _tol = 1e-4
    def __init__(self,model=None):
        """ Initialize. 

            Default uses the Maxwellian weighted Monte Carlo 
            sampler. Cross sections given in cm^2 by default.
            Standard energy range is 0-100 keV.

            Args:
                model: The physics model.
        """
 
        if model == None:
            self.model = MaxwellWeightedSampler( \
                           AstroModel(),InteractionModel())
            self.model.interaction.cross_section.total_xs = \
                1.0 * (units.cm)**2
            self.model.initialize()
        else:
            self.model = Model()

        self._Emin = 0
        self._Emax = 100 * units.keV



    def set_limits(self,Emin,Emax):
        """ Set the energy range the analysis looks at.
  
            Args:
                Emin: Minimum energy.
                Emax: Maximum energy.
        """
        self._Emin = Emin
        self._Emax = Emax

    def set_tolerance(self,tol):
        """ Set the tolerance used in weight and cross section
            calculations. Fractional value. 
      
            Args:
                tol: The new value
        """
        self._tol = tol

    def get_rate(self):
        """ Get the interaction rate.

            Here, we only need the overall rate the detector
            is expected to see for the given model.
        """
        err = 2
        rate = 0
        rate_var = 0
        N = 0 
        while err > _tol:
            s = model.sample()
            if self._Emin <= Er < self._Emax and s.weight>0:
                rate = rate + self.weight
                rate_var = rate_var + self.weight*self.weight
                err = np.sqrt(rate_var) / rate
            N = N + 1
        rate = rate / N
        return rate

    def get_limit(self,Mx,N_exp = 0, CL=0.9):
        """ Get the CL upper limit
            Default is 0 events with a 90% confidence level.

            Args:
                Mx: The WIMP mass
                N_exp: (int) The number of events measured
                CL: (float, 0-1) The confidence level.
        """

        self.model.interaction.Mx = Mx
        rate = self.get_rate()

        n_limit, fn = poisson_freq_ul(N_exp, CL)
        xs_limit = n_limit * \
                   self.model.cross_section.total_xs / rate
        return xs_limit, fn

def upper_limit(N_exp=0,CL=0.9):
    """ Function to calculate Poisson upper limits.

        Args:
            N_exp: Number of measured events.
            CL: Confidence level
    """
    # We measured nothing. Easy case
    if N_exp == 0:
        # Works for frequentist and also 
        # Bayesian with flat prior
        # The limit in number of events
        n_limit = -np.log(1 - CL) 
        return n_limit

    # We actually measured some events
    err = 2
    fn_last = np.zeros(3) + 2
    n_limit = 0.5*np.sqrt(N_exp)+N_exp # First guess
    fn = 0
    while err > UpperLimitBkgFree._tol:
        fn = ( CL - 1)
        dfdn = 0
        for i in range(N_exp+1): # Should the +1 be here?
            pmf = scipy.stats.poisson.pmf(i,n_limit)
            fn = fn + pmf 
            dfdn = dfdn + (-1 + i / n_limit) * pmf
        ## Newton's method
        

        #dfdn = np.sign(dfdn) * np.max([np.abs(dfdn),0.01])
        n_limit = n_limit - fn / dfdn
        if n_limit <= 0:
            print('Error: Mean has fallen below 0')
            print(n_limit,fn,dfdn,fn/dfdn)
            return 0
        err = np.max( np.abs(fn )) / (1-CL)
        if np.abs( 1 - fn_last[2]/fn ) < (UpperLimitBkgFree._tol)**2 :
            n_limit = n_limit + UpperLimitBkgFree._tol
        fn_last[0] = fn_last[1]
        fn_last[1] = fn_last[2]
        fn_last[2] = fn
        
    return 0,n_limit
