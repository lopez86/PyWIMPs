from . import VelocityDist
from .. import units
import numpy as np
class AstroModel:

    def __init__(self,norm=False):
        self.set_random(np.random)
        self.velocity = VelocityDist()
        self.wimp_density = 0.3 * units.GeV / (units.cm**3)
        self.vE = 220 * units.km / units.sec * np.array([0,0,1])
        self.v0 = 220 * units.km / units.sec
        self.vesc = 550 * units.km / units.sec 
        if norm:
            self.velocity.normalize()

    @property
    def random(self):
        return self._rand

    @random.setter
    def set_random(self,r):
        self.rand = r
        self.velocity.set_random(r)

    def set_params(self,pars):
        if 'vE' in pars.keys():
            self.vE = pars['vE']
        if 'v0' in pars.keys():
            self.v0 = pars['v0']
        if 'vesc' in pars.keys():
            self.vesc = pars['vesc']
        if 'rhox' in pars.keys():
            self.wimp_density = pars['rhox']
        self.velocity.set_params(pars)
        
    @property
    def vE(self):
        return self._vE
 
    @vE.setter
    def set_vE(self,ve):
        self._vE = ve
        self.velocity.vE = ve

    @property
    def v0(self):
        return self._v0

    @v0.setter
    def set_v0(self,v):
        self._v0 = v
        self.velocity.v0 = v

    @property
    def vesc(self):
        return self._vesc

    @vesc.setter
    def set_vesc(self,v):
        self._vesc = v
        self.velocity.vesc = v

    @property
    def wimp_density(self):
        return self._wimp_density

    @wimp_density.setter
    def set_wimp_density(self,rho):
        self._wimp_density = rho

    @property
    def rho(self):
        return self._wimp_density

