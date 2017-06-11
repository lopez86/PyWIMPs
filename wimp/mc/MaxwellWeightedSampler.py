import numpy as np
from .sample import Sample
from .. import mathtools
from .. import units
class MaxwellWeightedSampler:
    def __init__(self,astro_model,int_model):
        self.astro_model = astro_model
        self.interaction = int_model

    def set_params(self,pars,set_models=False):
        if set_models:
            self.astro_model.set_params(pars)
            self.interaction.set_params(pars)

    def initialize(self):
        self.vE = self.astro_model.vE()
        self.v0 = self.astro_model.v0()
        self.vesc = self.astro_model.vesc()
        self.Mx = self.interaction.Mx()
        self.Mt = self.interaction.Mt()
        self.xs = self.interaction.total_xs()
        self.mu = self.Mx*self.Mt / (self.Mx+self.Mt)
        self.Mtot = self.interaction.Mtot()
        self.rho = self.astro_model.wimp_density()
        self.e1,self.e2,self.e3 = mathtools.get_axes(self.vE)
        self.vE_mag = np.sqrt(self.vE.dot(self.vE))

    def sample(self):


        vec = np.random.normal(-self.vE,self.v0/np.sqrt(2),3)
        vec2 = vec + self.vE

            #Probability to throw this from the full Maxwellian distribution
        vec_prob = 1./(np.pi*self.v0*self.v0)**1.5 * \
                   np.exp( - (vec2.dot(vec2)) / (self.v0*self.v0) ) 

        vec_mag = np.sqrt(vec.dot(vec))
        # The factor of 1./vec_mag comes from v from the flux and 1/v^2
        # from the recoil energy normalization
        Ex = 0.5 * self.Mx * (vec_mag/units.speed_of_light)**2
        Emax = self.interaction.cross_section.MaxEr(Ex)
        weight = self.astro_model.velocity.f(vec) / vec_prob * vec_mag 
        #return Sample(0,0,weight,0)
        ## Now let's look at the interaction part

        E = np.random.rand() * Emax
        phi = np.random.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosThetaLab(Ex,E)

        ## Add in the form factor
        Q2 = 2 * self.Mt * E
        weight = weight * self.interaction.form_factor.ff2(Q2)

        # Add in the constants so that we normalize to rate
       # weight = weight * self.xs * self.rho * self.Mtot / (2 * self.mu*self.mu * self.Mx) * Emax
        weight *= self.xs * self.rho/self.Mx * self.Mtot/self.Mt
        
        ## Let's go back into the lab frame:
        e1v,e2v,e3v = mathtools.get_axes(vec)
        recoil_lab = e1v * cosTheta + \
                     np.sqrt(1-cosTheta*cosTheta) * \
                     ( np.cos(phi) * e2v + np.sin(phi) * e3v )

        return Sample(E,recoil_lab,weight,vec)

