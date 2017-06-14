""" UniformWeightedSampler.py

    Monte Carlo sampling drawing from a uniform
    distribution to get samples with generator-level 
    weights.
"""
__author__ = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np
from .. import units
from .. import mathtools
from .sample import Sample


class UniformWeightedSampler:
    """ Class to perform sampling on the standard
        halo and cross section model. Throws are based on
        a uniform distribution, so biased samples
        are obtained along with appropriate weights to get
        the correct distribution.

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
            e1min: Minimum velocity in the e1 direction
            e1max: Maximum velocity in the e1 direction
            e2min: Minimum velocity in the e2 direction
            e2max: Maximum velocity in the e2 direction
            e3min: Minimum velocity in the e3 direction
            e3max: Maximum velocity in the e3 direction
    """
    def __init__(self,astro_model,int_model):
        """ Initialize the object. 

            Args:
                astro_model (AstroModel)
                int_model (InteractionModel)
        """
        self.astro_model = astro_model
        self.interaction = int_model
        self._rand = np.random

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
        """
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
        self.xs = self.interaction.total_xs
        self.mu = self.Mx*self.Mt / (self.Mx+self.Mt)
        self.Mtot = self.interaction.Mtot
        self.rho = self.astro_model.wimp_density

        self.e1,self.e2,self.e3 = mathtools.get_axes(self.vE)

        self.vE_mag = np.sqrt(self.vE.dot(self.vE))
        self.e1min = -self.vesc - self.vE_mag
        self.e1max = self.vesc - self.vE_mag
        self.e2min = -self.vesc
        self.e2max = self.vesc
        self.e3min = -self.vesc
        self.e3max = self.vesc

        self.vol = ((self.e1max-self.e1min)
                   * (self.e2max-self.e2min) 
                   * (self.e3max-self.e3min))


    def sample(self):
        """ Get a sample. 

            Returns:
                A biased sample with a generator weight
        """

        rnd = self._rand.rand(3)
        vec = self.e1*((self.e1max-self.e1min) * rnd[0] +self.e1min) 
        vec = vec + self.e2*((self.e2max-self.e2min) * rnd[1] + self.e2min)
        vec = vec + self.e3*((self.e3max-self.e3min) * rnd[2] + self.e3min)
        vec_mag = np.sqrt(vec.dot(vec))
        # The factor of 1./vec_mag comes from v from the flux and 1/v^2
        # from the recoil energy normalization
        weight = self.vol * self.astro_model.velocity.f(vec) * vec_mag
        
        ## Now let's look at the interaction part
        Ex = 0.5 * self.Mx * (vec_mag/units.speed_of_light)**2
        Emax = self.interaction.cross_section.MaxEr(Ex)

        E = self._rand.rand() * Emax
        phi = self._rand.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosThetaLab(Ex,E)
        ## E and phi are uniform so they do not introduce a weight
        ## That is, normalization and volume cancel out

        ## Add in the form factor
        Q2 = 2 * self.Mt * E
        weight = weight * self.interaction.form_factor.ff2(Q2)
        # Add in the constants so that we normalize to rate
        
        weight *= self.xs * self.rho/self.Mx * self.Mtot/self.Mt 


        ## Let's go back into the lab frame:
        e1v,e2v,e3v = mathtools.get_axes(vec)
        recoil_lab = (e1v * cosTheta 
                      + np.sqrt(1-cosTheta*cosTheta) 
                      * ( np.cos(phi) * e2v + np.sin(phi) * e3v ))

        return Sample(E,recoil_lab,weight,vec)

