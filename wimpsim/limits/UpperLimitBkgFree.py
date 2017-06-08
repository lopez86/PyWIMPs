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

    _tol = 1e-4
    def __init__(self,model=None):
 
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
        self._Emin = Emin
        self._Emax = Emax

    def set_tolerance(self,tol):
        self._tol = tol

    def get_rate(self):
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

        self.model.interaction.Mx = Mx
        rate = self.get_rate()

        n_limit, fn = upper_limit(N_exp, CL)
        xs_limit = n_limit * \
                   self.model.cross_section.total_xs / rate
        return xs_limit, fn

def upper_limit(N_exp=0,CL=0.9):
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
        
    return n_limit
