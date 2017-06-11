from . import VelocityDist
from .. import units
import numpy as np
class AstroModel:

    def __init__(self,norm=False):
        self.rand = np.random
        self.velocity = VelocityDist()
        self._wimp_density = 0.3 * units.GeV / (units.cm**3)
        self._vE = 220 * units.km / units.sec * np.array([0,0,1])
        self._v0 = 220 * units.km / units.sec
        self._vesc = 550 * units.km / units.sec 
        self.fill_params()
        if norm:
            self.velocity.normalize()        
    def set_random(self,r):
        self.rand = r
        self.velocity.set_random(r)

    def fill_params(self):
        pars = {'rhox':self._wimp_density,'vE':self._vE,'v0':self._v0,'vesc':self._vesc}
        self.velocity.set_params(pars)

    def set_params(self,pars):
        if 'vE' in pars.keys():
            self._vE = pars['vE']
        if 'v0' in pars.keys():
            self._v0 = pars['v0']
        if 'vesc' in pars.keys():
            self._vesc = pars['vesc']
        if 'rhox' in pars.keys():
            self._wimp_density = pars['rhox']
        self.velocity.set_params(pars)
        

    def vE(self):
        return self._vE
 
    def v0(self):
        return self._v0

    def vesc(self):
        return self._vesc

    def wimp_density(self):
        return self._wimp_density

