import numpy as np
from .sample import Sample
from .. import mathtools
from .. import units
class AcceptRejectSampler:
    def __init__(self,astro_model,int_model):
        self.rand = np.random
        self.astro_model = astro_model
        self.interaction = int_model
        self.max_iter = 10000

    def set_random(self,r,set_models=False):
        self.rand = r
        if set_models:
            self.astro_model.set_random(r)
            self.interaction.set_random(r)

    def set_params(self,pars,set_models=False):
        if 'AccRejMaxIter' in pars.keys():
            pars.max_iter = pars['AccRejMaxIter']
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

        #Velocity that maximizes the Maxwell-Boltzmann distribution
        #Ignores any effects from the truncation
        #self.vmaxP =  0.5 * (self.vE_mag + np.sqrt(self.vE_mag*self.vE_mag + 2*self.v0*self.v0))
        self.vmaxP =  0.5 * (self.vE_mag + np.sqrt(self.vE_mag*self.vE_mag + 6*self.v0*self.v0))

        if self.vE_mag<1e-12:
            self.vmaxP_vec = -self.vmaxP * np.array([0,0,1])
        else:
            self.vmaxP_vec = -self.vmaxP * self.vE / self.vE_mag
        #self.maxP = self.vmaxP * self.astro_model.velocity.f_no_escape(self.vmaxP_vec)
        self.maxP = self.vmaxP**3 * self.astro_model.velocity.f_no_escape(self.vmaxP_vec)
        # Maximum possible WIMP velocity
        self.vmax = self.vesc + self.vE_mag
        # Minimum possible WIMP velocity
        self.vmin = -min(self.vesc - self.vE_mag,0)


    def sample(self):


        passed = False

        vec = np.array([0,0,0])
        Er = -1       
        iteration = 0
        while passed is False:

            if (iteration >= self.max_iter):
                print("Accept/Reject: Max iteration reached")
                return(-1,np.array([0,0,0]),0)
            # First, throw a velocity:
            v = self.rand.rand() * (self.vmax-self.vmin) + self.vmin
            cosTh = 2 * self.rand.rand() - 1
            phi = self.rand.rand() * 2*np.pi
            sinTh = np.sqrt(1-cosTh*cosTh)
            vec = np.array([v*sinTh*np.cos(phi),v*sinTh*np.sin(phi),v*cosTh])
            vec_mag = np.sqrt(vec.dot(vec))
            Ex = 0.5 * self.Mx * (vec_mag/units.speed_of_light)**2
            Emax = self.interaction.cross_section.MaxEr(Ex)
            # Throw a recoil energy
            Er = self.rand.rand() * Emax
            Q2 = 2 * self.Mt * Er

            # Throw a random probability
            rnd = self.rand.rand() * self.maxP

            # Calculate the probability:
            P = vec_mag**3 * self.astro_model.velocity.f(vec) * \
                self.interaction.form_factor.ff2(Q2)
            if P > self.maxP: 
              print( 'Illegal P found: ' +str(P/self.maxP))
            # Compare:
            if P > rnd:
              break

            iteration = iteration + 1

        phi = self.rand.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosThetaLab(Ex,Er)


        ## Let's go back into the lab frame:
        e1v,e2v,e3v = mathtools.get_axes(vec)
        recoil_lab = e1v * cosTheta + \
                     np.sqrt(1-cosTheta*cosTheta) * \
                     ( np.cos(phi) * e2v + np.sin(phi) * e3v )

        return Sample(Er,recoil_lab,1,vec)

