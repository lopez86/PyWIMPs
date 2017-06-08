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

    _tol = 1e-2
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

    def get_interval(self,Mx,N_exp = 0,N_bkg=0, CL=0.9):

        self.model.interaction.Mx = Mx
        rate = self.get_rate()

        n_limit, fn = fc_interval(N_exp,N_bkg, CL)
        xs_limit = n_limit * \
                   self.model.cross_section.total_xs / rate

        return xs_min,xs_max

def fc_limits(s,N_exp=0,b=0,CL=0.9):

    # Get approx. # of sigma:
    sigma = scipy.special.erfinv(CL) * np.sqrt(2)
    N_max = int(max(N_exp + 4*sigma * np.sqrt(N_exp),20) )
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
    
    low_limit = 0
    sigma = scipy.special.erfinv(CL) * np.sqrt(2)
    up_limit = int(max(N_exp+10*sigma*np.sqrt(N_exp),20))

    lim_minl,lim_maxl,total_probl = fc_limits(low_limit,N_exp,b,CL)

    lim_minu,lim_maxu,total_probu = fc_limits(up_limit,N_exp,b,CL)
    while lim_minu <= N_exp:
        up_limit = 2 * up_limit
        lim_minu,lim_maxu,total_probu = \
                 fc_limits(up_limit,N_exp,b,CL)

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
            minl, maxl, total_prob = fc_limits(new_l,N_exp,b,CL)
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
        minl, maxl, total_prob = fc_limits(new_l,N_exp,b,CL)
        if minl <= N_exp:
            low_l = new_l
        else:
            up_l = new_l
    up_limit = low_l


    return low_limit,up_limit
