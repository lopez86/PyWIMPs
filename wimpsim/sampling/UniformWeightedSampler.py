import numpy as np
from .. import units
from .sample import Sample

class UniformWeightedSampler:
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
        self.e1min = -self.vesc - self.vE_mag
        self.e1max = self.vesc - self.vE_mag
        self.e2min = -self.vesc
        self.e2max = self.vesc
        self.e2min = -self.vesc
        self.e2max = self.vesc

        self.vol = (self.e1max-self.e1min)*(self.e2max-self.e2min)*(self.e3max-self.e3min)


    def sample():


        rnd = np.random.rand(3)
        vec = (self.e1max-self.e1min) * rnd[0] +self.e1min
        vec = vec + (self.e2max-self.e2min) * rnd[1] + self.e2min
        vec = vec + (self.e3max-self.e3min) * rnd[2] + self.e3min
        vec_mag = sqrt(vec.dot(vec))
        # The factor of 1./vec_mag comes from v from the flux and 1/v^2
        # from the recoil energy normalization
        weight = self.vol * self.astro.velocity.f(self.vec) / self.vec_mag

        ## Now let's look at the interaction part
        Ex = 0.5 * self.Mx * (vec/units.speed_of_light)**2
        Emax = self.interaction.cross_section.Emax(Ex)

        E = np.random.rand() * Emax
        phi = np.random.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosTheta(E)
        ## E and phi are uniform so they do not introduce a weight
        ## That is, normalization and volume cancel out

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

