""" mcmc_tuning.py

Shows some diagnostics when trying to tune the Markov Chain
sampler.

Requires PyROOT

Example:

>>> run_test()

will test by default a burnin period of 5000 samples and 
a proposal distribution with a standard deviation (1-dim)
of 100 km/s.

"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from .. import units
from ..mc.MCMCSampler import MCMCSampler
from ..mc.MaxwellWeightedSampler import MaxwellWeightedSampler
from ..astro.AstroModel import AstroModel
from ..xsec.InteractionModel import InteractionModel
from ROOT import TH1F

import numpy as np
import matplotlib.pyplot as plt


def get_test_hist(am,im):
    """ Get a histogram to use to test the MCMC result.
   
    Uses 4M throws from the Maxwell-Boltzmann weighted 
    sampler. 

    Args:
        am: The astrophysics model
        im: The interaction model

    Returns:
        A ROOT histogram showing the energy spectrum

    """
    maxw = MaxwellWeightedSampler(am,im)
    maxw.initialize()
    h = TH1F('maxw','Maxwell;E_{r} [keV];Rate [shape only]',60,0,120)
    h.Sumw2()
    n = 4000000
    for i in range(n):
        if i % (n//10) == 0:
            print('Test histogram %0.1f%% done' %( (10.0*i) / (n//10) ) )
        ms = maxw.sample()
        h.Fill(ms.Er / units.keV,ms.weight)
    h.Scale(1./n)
    return h


def evaluate_params(am,im,nburnin,sigma,hmaxw, N,ns):
    """ Evaluate how the starting parameters perform
        when comparing results to the reference histogram.

    Args:
        am: The astrophysics model
        im: The interaction model
        nburnin: The length of the burn-in period
        sigma: The width of the proposal distribution
        hmaxw: The reference histogram
        N: The number of samples to generate
        ns: The number of samples to average for diagnostics
    
    Returns:
        samples: The sample points at which we get diagnostics.
                 These are averaged over ns individual samples
        accept_frac: The acceptance rate of throws around each
                     sample point. Should be stable if things
                     are working well.
        chisq: The chisq difference between the MCMC histogram
               and the reference histogram. Will likely get 
               worse with huge numbers of events as the uncertainties
               of the MCMC get better than the weighted throws.
               The errors here may still need some work for this
               to work well.
    """
    mcmc = MCMCSampler(am,im)
    pars={'MCMCburnin':nburnin,'MCMCSigma':sigma }

    mcmc.set_params(pars)
    mcmc.initialize()

    ntries = 0
    #hmcmc = TH1F('MCMC','Metropolis;E_{r} [keV];Rate [shape only]',60,0,120)
    hmcmc = hmaxw.Clone('MCMC')
    hmcmc.SetTitle('Metropolis')
    
    accept_frac = np.zeros(N//ns)
    chisq = np.zeros(N//ns)
    samples = ns * np.array(range(N//ns))

    for i in range(N):
        if i % (N//10) == 0:
            print('MCMC %0.1f%% done' %( (10.0*i) / (N//10) ) )
        
        if (i % ns) == (ns-1):
            accept_frac[i//ns] = (1.0*ns) / ntries        
            ntries = 0
            chisquare = 0
            sum_ratio = hmaxw.Integral() / hmcmc.Integral()
            for j in range(1,hmcmc.GetNbinsX()+1):
                nev = hmcmc.GetBinContent(j)
                estimate = hmaxw.GetBinContent(j) / sum_ratio
                esterr = hmaxw.GetBinError(j) / sum_ratio
                # Ignore the uncertainty on hmaxw

                chisq[i//ns] = chisq[i//ns] +  \
                               (nev - estimate )**2 / (estimate + esterr**2)
            chisq[i//ns] = chisq[i//ns] / (hmaxw.GetNbinsX()) 
            # Approx. Chi^2/NDF

        sample = mcmc.sample()
        hmcmc.Fill(sample.Er/units.keV)
        ntries = ntries + mcmc.Ntries

    
    return samples, accept_frac, chisq

def run_test(nburnin_vec=[5000,],sigma_vec=[100*units.km/units.sec,]):
    """ Loop over a bunch of different metaparameter values.

        Currently, only the last one is used for anything, so length 1
        lists are probably best. After done looping, makes some plots
        about the last test done.

        Args:
            nburnin_vec: List of burn-in lengths to check
            sigma_vec: List of proposal density widths to check
    """
    am = AstroModel()
    im = InteractionModel()

    hmaxw = get_test_hist(am,im)
    for nb in nburnin_vec:
        for sig in sigma_vec:
            samples, accept_frac, chisq =  \
                     evaluate_params(am,im,nb,sig,hmaxw,500000,500)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    ax.plot(samples,chisq)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Histogram ChiSquare')
#    ax.set_ylim(0,10)
    fig2 = plt.figure(2)
    ax = fig2.add_subplot(111)
    ax.plot(samples,accept_frac)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Acceptance Fraction')

    plt.show()

