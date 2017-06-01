import numpy as np
from .sample import Sample

class MaxwellWeightedSampler:
    def __init__(self,astro_model,int_model):
        self.astro = astro_model
        self.interaction = int_model

    def initialize():
        self.vE = self.astro.vE()
        self.v0 = self.astro.v0()
        self.vesc = astro.vesc()
        self.Mx = self.interaction.Mx()
        self.Mt = self.interaction.Mt()
        self.xs = self.interaction.total_xs()
        self.mu = self.Mx*self.Mt / (self.Mx+self.Mt)
        self.Mtot = self.interaction.Mtot()
        self.rho = self.astro.wimp_density()
        self.e1,self.e2,self.e3 = wsmath.get_axes(self.vE)
        self.vE_mag = np.sqrt(self.vE.dot(self.vE))

    def sample():


        vec = np.random.normal(-self.vE,self.v0/np.sqrt(2),3)
        vec2 = vec + self.vE

            #Probability to throw this from the full Maxwellian distribution
        vec_prob = 1./(np.pi*self.v0*self.v0)**1.5 * \
                   np.exp( - (vec2.dot(vec2)) / (self.v0*self.v0) ) 

        vec_mag = np.sqrt(vec.dot(vec))
        # The factor of 1./vec_mag comes from v from the flux and 1/v^2
        # from the recoil energy normalization
        weight = self.astro.velocity.f(vec) / vec_mag / vec_prob

        ## Now let's look at the interaction part
        Ex = 0.5 * self.Mx * (vec/units.speed_of_light)**2
        Emax = self.interaction.cross_section.Emax(Ex)

        E = np.random.rand() * Emax
        phi = np.random.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosTheta(E)

        ## Add in the form factor
        Q2 = 2 * self.Mt * E
        weight = weight * self.interaction.form_factor.ff2(Q2)

        # Add in the constants so that we normalize to rate
        weight = weight * self.xs * self.rho * self.Mtot / (2 * self.mu*self.mu * self.Mx) * Emax

        ## Let's go back into the lab frame:
        e1v,e2v,e3v = wsmath.get_axes(vec)
        recoil_lab = e1v * cosTheta + \
                     np.sqrt(1-cosTheta*cosTheta) * \
                     ( np.cos(phi) * e2v + np.sin(phi) * e3v )

        return Sample(E,recoil_lab,weight,vec)

