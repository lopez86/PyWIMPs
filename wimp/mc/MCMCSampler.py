""" MaxwellWeightedSampler.py

    Monte Carlo sampling with a Markov Chain using the
    Metropolis-Hastings algorithm to get samples that are
    unbiased when large datasets are generated.
"""
__author__ = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np
from .. import units
from .. import mathtools
from .sample import Sample
import numpy as np


class MCMCSampler:
    """ Class to perform sampling on the standard
        halo and cross section model. Throws are based on
        a Markov Chain Monte Carlo using the Metropolis-Hastings 
        algorithm with a Gaussian proposal distribution.

        For large numbers of samples, the samples are unbiased.
 
        The user should tune the sampler to make sure that the
        results are working properly for the parameters chosen
        and the number of samples being generated.

        If you are not sure how to evaluate the performance,
        use another sampler such as AcceptRejectSampler instead.

        Attributes:
            astro_model (AstroModel)
            interaction (InteractionModel)
            max_iter (int): Max # of iterations before
                            stopping throws

            vE: Earth velocity vector
            v0: Dispersion velocity
            vesc: Galactic escape velocity
            Mx: Dark matter mass
            Mt: Target nucleus mass
            xs: WIMP-nucleus cross section
            mu: Interaction reduced mass
            Mtot: Total detector mass
            rho: WIMP mass density
            e1: Unit vector along vE
            e2: Unit vector orthogonal to vE
            e3: Unit vector orthogonal to vE
            vE_mag: vE magnitude
            sigma: The width (1D) of the Gaussian used
                   to throw velocities.
            nburnin: The number of samples to get
                     on initialization
            Ntries: The number of throws used to get
                    the last sample.
            lastv: The WIMP velocity for the last sample
            lastE: The recoil energy for the last sample
            lastP: The probability density for the last sample
            Ex: The WIMP energy for the current sample
 
    """
    def __init__(self,astro_model,int_model):
        self._rand = np.random
        self.astro_model = astro_model
        self.interaction = int_model
        self.sigma = 20 * units.km / units.sec
        self.nburnin = 1000
        self.Ntries = 0

    @property
    def random(self):
        """ Random number generator. """
        return self._rand

    @random.setter
    def set_random(self,r,set_models=False):
        """ Set the random number generator.
     
            Args:
                r (Numpy RandomState)
                set_models: Also set the generator for
                            the models
        """
        self._rand = r
        if set_models:
            self.astro_model.set_random(r)
            self.interaction.set_random(r)

    def set_params(self,pars,set_models=False):
        """Set the parameters based on a dictionary.

           Args:
               pars {string}
               set_models: Also set the parameters for
                           the models

           Parameters:
               MCMCSigma: The 1D width of the proposal velocity
                          distribution
               MCMCNburnin: The number of throws to use on 
                            initialization to avoid bias
        """
        if 'MCMCSigma' in pars.keys():
            self.sigma = pars['MCMCSigma']
        if 'MCMCNburnin' in pars.keys():
            self.nburnin = pars['MCMCnburnin']
        if set_models:
            self.astro_model.set_params(pars)
            self.interaction.set_params(pars)

    def initialize(self):
        """ Perform a final initialization to prepare for
            generating samples.
        """
        self.vE = self.astro_model.vE
        self.v0 = self.astro_model.v0
        self.vesc = self.astro_model.vesc
        self.Mx = self.interaction.Mx
        self.Mt = self.interaction.Mt
        self.mu = self.Mx*self.Mt / (self.Mx+self.Mt)

        self.e1,self.e2,self.e3 = mathtools.get_axes(self.vE)

        self.vE_mag = np.sqrt(self.vE.dot(self.vE))

        ### Now, run the burnin part

        ### First, set up an initial guess:
        vguess = 0.95*self.vesc * (self._rand.rand())**(0.33333)        
        costhguess = 2 * self._rand.rand() - 1
        phiguess = 2*np.pi * self._rand.rand()
        sinthguess = np.sqrt(1 - costhguess*costhguess)
        self.lastv = np.array([vguess * np.cos(phiguess) * sinthguess,
                               vguess * np.sin(phiguess) * sinthguess,
                               vguess * costhguess] )
        self.lastv = self.lastv - self.vE
        vguess = np.sqrt(self.lastv.dot(self.lastv))
        Ex = 0.5 * self.Mx * (vguess/units.speed_of_light)**2
        Emax = self.interaction.cross_section.MaxEr(Ex)
        self.lastE = self._rand.rand() * Emax
        Q2 = 2 * self.Mt * self.lastE        

        self.lastP = (vguess
                      * self.astro_model.velocity.f(self.lastv)
                      * self.interaction.form_factor.ff2(Q2))
                     
        # The other things that we're throwing (Er, phi_r)
        # will have flat priors so we don't need to save anything

        for idx in range(self.nburnin): 
            self.sample()

 
    def sample(self):
        """ Get a sample. Be careful: For a Markov Chain
            nearby samples are correlated. For uncorrelated throws,
            you must get a dataset with many samples and then draw 
            randomly from that dataset.
            

            Returns:
                An unbiased sample
        """

        self.Ntries = 0
        done = False
        Ex = 0
        while not done:
            self.Ntries = self.Ntries + 1
            # Propose a new WIMP velocity
            vprop = self._rand.normal(self.lastv,self.sigma)
            vec_mag = np.sqrt(vprop.dot(vprop))
            ## Get the maximum energy
            Ex = 0.5 * self.Mx * (vec_mag/units.speed_of_light)**2
            Emax = self.interaction.cross_section.MaxEr(Ex)

            ## There is a 1/Emax from the differential cross section
            ## But:

            ## MH: q(last | this) p(this) /   q(this | last ) p(last)
            ## The velocity throws are symmetric and phi doesn't matter here
            ## 
            ## --> q(Elast | Ethis) p(this) / q(Ethis | Elast) p(last)
            ##  --> Emaxthis / Emaxlast *  p(this) / p(last)
            ## p ~ vf(v)|F|^2 / Emax
            ## --> vf(v)|F|^2
            ## Acceptance only depends on E through the form factor

            ## Propose an energy:
            Eprop = self._rand.rand()*Emax
            Q2 = 2 * self.Mt * Eprop
 
            ## Get the probability
            Pprop = (vec_mag
                         * self.astro_model.velocity.f(vprop)
                         * self.interaction.form_factor.ff2(Q2))

           
            if Pprop <= 0:
                continue # Don't need to waste time
            # Get the acceptance function:
            alpha = min(1,Pprop / self.lastP)
            # Throw a random number
            if alpha <= self._rand.rand():
                continue

            self.lastv = vprop
            self.lastE = Eprop
            self.lastP = Pprop
            break

        phi = self._rand.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosThetaLab(Ex,
                            self.lastE)

        ## Let's go back into the lab frame:
        e1v,e2v,e3v = mathtools.get_axes(self.lastv)
        recoil_lab = (e1v * cosTheta
                      + np.sqrt(1-cosTheta*cosTheta)
                      * ( np.cos(phi) * e2v + np.sin(phi) * e3v ))

        return Sample(self.lastE,recoil_lab,1,self.lastv)

