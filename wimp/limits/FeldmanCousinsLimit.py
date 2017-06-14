""" FeldmanCousinsLimit.py

Calculate limits and confidence intervals using the 
Feldman-Cousins method to get proper coverage.

Used to get confidence intervals for signals given
some background and an experiment typically with
a small number of selected events.
"""
__author__ =    'Jeremy P. Lopez'
__date__ =      'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'


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
import scipy.special


class FeldmanCousinsLimit:
    """ Class to calculate Feldman-Cousins confidence intervals.
 
        Feldman and Cousins published a method to get proper 
        coverage when switching from limits to bounded intervals
        by introducing an ordering scheme when creating an
        interval.

        Attributes:
            _tol: Default: 1%. Maximum bound on the estimated
                  error of rate calculations. Also used to get
                  the maximum uncertainty on the interval limits.
                  static class member.

            model: Model used to calculate rates. Must have
                   sample() defined.
            _Emin: Minimum energy
            _Emax: Maximum energy 

    """
    _tol = 1e-2
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

    def get_interval(self,Mx,N_exp = 0,N_bkg=0, CL=0.9):
        """ Get the confidence interval.
            Default is 0 events with a 90% confidence level.

            Args:
                Mx: The WIMP mass
                N_exp: (int) The number of events measured
                N_bkg: (float) The estimated number of backgrounds
                CL: (float, 0-1) The confidence level.
        """
        self.model.interaction.Mx = Mx
        rate = self.get_rate()

        n_limit, fn = fc_interval(N_exp,N_bkg, CL)
        xs_limit = n_limit * \
                   self.model.cross_section.total_xs / rate

        return xs_min,xs_max

def fc_limits(s,b=0,CL=0.9):
    """ Function to calculate bounds on the number of measured
        events for a model given the Feldman-Cousins ordering
        principle.

        Args:
            s: The expected number of signal events
            b: The expected number of background events
            CL: The confidence interval

        Returns:
            The CL bounds on the number of expected events
            in the experiment given a signal and background
            model.
    """
    # Get approx. # of sigma:
    sigma = scipy.special.erfinv(CL) * np.sqrt(2)
    N_max = int(max(s+b + 4*sigma * np.sqrt(s+b),20) )
    # Should be enough for a reasonable data set
    while True:
        N_max = N_max * 2
        ratio = np.zeros(N_max)
        pval = np.zeros(N_max)
        for n in range(N_max):
            pval[n] = scipy.stats.poisson.pmf(n,s+b)
            pbest = 1e4
            if b > n:
                pbest = scipy.stats.poisson.pmf(n,b)
          
                # minus sign is because sorting is in ascending order
            else:
                pbest = scipy.stats.poisson.pmf(n,n)
            ratio[n] = -pval[n] / pbest  
        sorted_indices = np.argsort(ratio)
        total_prob = 0
        lim_min = N_max
        lim_max = -1
        n = 0
        while total_prob < CL and n < len(ratio):
            total_prob += pval[sorted_indices[n]]
            if sorted_indices[n] < lim_min:
                lim_min = sorted_indices[n]
            if sorted_indices[n] > lim_max:
                lim_max = sorted_indices[n]
            n = n+1
    
        if lim_max < N_max-1:
            break
  
    return lim_min,lim_max,total_prob

def fc_interval(N_exp=0,b=0,CL=0.9):
    """ Function to calculate Feldman-Cousins confidence
        intervals.

        Args:
            N_exp: (int)  Number of measured events. Default: 0
            b: (float) Expected number of backgrounds. Default: 0
            CL: (float, 0-1) Confidence level. 
    """
    low_limit = 0
    sigma = scipy.special.erfinv(CL) * np.sqrt(2)
    up_limit = int(max(N_exp+10*sigma*np.sqrt(N_exp),20))

    lim_minl,lim_maxl,total_probl = fc_limits(low_limit,b,CL)

    lim_minu,lim_maxu,total_probu = fc_limits(up_limit,b,CL)
    while lim_minu <= N_exp:
        up_limit = 2 * up_limit
        lim_minu,lim_maxu,total_probu = \
                 fc_limits(up_limit,b,CL)

    upper_lim = False
    
    if lim_minl <= N_exp < lim_maxl:
        upper_lim = True # We want an interval in this case
        print("Calculating upper limit")
    # Binary search for the limit
    if upper_lim == False:
        low_l = low_limit
        up_l = up_limit
        
        while up_l - low_l > FeldmanCousinsLimit._tol:
            new_l = 0.5 * (up_l - low_l) + low_l
            minl, maxl, total_prob = fc_limits(new_l,b,CL)
            if maxl >= N_exp:
                up_l = new_l
            else:
                low_l = new_l
        low_limit = low_l

    # Upper limit now

    low_l = 0
    up_l = up_limit
    
    while up_l - low_l > FeldmanCousinsLimit._tol:
        new_l = 0.5 * (up_l - low_l) + low_l
        minl, maxl, total_prob = fc_limits(new_l,Nb,CL)
        if minl <= N_exp:
            low_l = new_l
        else:
            up_l = new_l
    up_limit = low_l


    return low_limit,up_limit
