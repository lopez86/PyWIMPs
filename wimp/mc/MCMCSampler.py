import numpy as np
from .. import units
from .. import mathtools
from .sample import Sample


class MCMCSampler:
    def __init__(self,astro_model,int_model):
        self.astro_model = astro_model
        self.interaction = int_model
        self.sigma = 20 * units.km / units.sec
        self.nburnin = 1000
        self.Ntries = 0
    def set_params(self,pars,set_models=False):
        if 'MCMCSigma' in pars.keys():
            self.sigma = pars['MCMCSigma']
        if 'MCMCNburnin' in pars.keys():
            self.nburnin = pars['MCMCnburnin']
        if set_models:
            self.astro_model.set_params(pars)
            self.interaction.set_params(pars)

    def initialize(self):
        self.vE = self.astro_model.vE()
        self.v0 = self.astro_model.v0()
        self.vesc = self.astro_model.vesc()
        self.Mx = self.interaction.Mx()
        self.Mt = self.interaction.Mt()
        self.mu = self.Mx*self.Mt / (self.Mx+self.Mt)

        self.e1,self.e2,self.e3 = mathtools.get_axes(self.vE)

        self.vE_mag = np.sqrt(self.vE.dot(self.vE))
        self.e1min = -self.vesc - self.vE_mag
        self.e1max = self.vesc - self.vE_mag
        self.e2min = -self.vesc
        self.e2max = self.vesc
        self.e3min = -self.vesc
        self.e3max = self.vesc

        ### Now, run the burnin part

        ### First, set up an initial guess:
        vguess = 0.95*self.vesc * (np.random.rand())**(0.33333)        
        costhguess = 2 * np.random.rand() - 1
        phiguess = 2*np.pi * np.random.rand()
        sinthguess = np.sqrt(1 - costhguess*costhguess)
        self.lastv = np.array([vguess * np.cos(phiguess) * sinthguess, \
                               vguess * np.sin(phiguess) * sinthguess, \
                               vguess * costhguess] )
        self.lastv = self.lastv - self.vE
        vguess = np.sqrt(self.lastv.dot(self.lastv))
        Ex = 0.5 * self.Mx * (vguess/units.speed_of_light)**2
        Emax = self.interaction.cross_section.MaxEr(Ex)
        self.lastE = np.random.rand() * Emax
        Q2 = 2 * self.Mt * self.lastE        

        self.lastP = vguess * \
                     self.astro_model.velocity.f(self.lastv) * \
                     self.interaction.form_factor.ff2(Q2)
                     
        # The other things that we're throwing (Er, phi_r)
        # will have flat priors so we don't need to save anything

        for idx in range(self.nburnin): 
            self.sample()

 
    def sample(self):

        self.Ntries = 0
        done = False
        while not done:
            self.Ntries = self.Ntries + 1
            # Propose a new WIMP velocity
            vprop = np.random.normal(self.lastv,self.sigma)
            vec_mag = np.sqrt(vprop.dot(vprop))
            ## Get the maximum energy
            self.Ex = 0.5 * self.Mx * (vec_mag/units.speed_of_light)**2
            Emax = self.interaction.cross_section.MaxEr(self.Ex)

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
            Eprop = np.random.rand()*Emax
            Q2 = 2 * self.Mt * Eprop
 
            ## Get the probability
            Pprop = vec_mag * \
                         self.astro_model.velocity.f(vprop) * \
                         self.interaction.form_factor.ff2(Q2)

           
            if Pprop <= 0:
                continue # Don't need to waste time
            # Get the acceptance function:
            alpha = min(1,Pprop / self.lastP)
            # Throw a random number
            if alpha <= np.random.rand():
                continue

            self.lastv = vprop
            self.lastE = Eprop
            self.lastP = Pprop
            break

        phi = np.random.rand() * 2 * np.pi
        cosTheta = self.interaction.cross_section.cosThetaLab(self.Ex,self.lastE)

        ## Let's go back into the lab frame:
        e1v,e2v,e3v = mathtools.get_axes(self.lastv)
        recoil_lab = e1v * cosTheta + \
                     np.sqrt(1-cosTheta*cosTheta) * \
                     ( np.cos(phi) * e2v + np.sin(phi) * e3v )

        return Sample(self.lastE,recoil_lab,1,self.lastv)

