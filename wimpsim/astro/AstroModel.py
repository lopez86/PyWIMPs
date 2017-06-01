import wimpsim.astro.VelocityDist
import wimpsim.Units as units

class AstroModel:

    def __init__(self,rho):
        self._velocity = wimpsim.astro.VelocityDist()
        self._wimp_density = 0.3 * units.GeV / (units.cm**2)
        self._vE = 220 * units.km / units.sec
        self._v0 = 230 * units.km / units.sec
        self._vesc = 550 * units.km / units.sec 
        self.fill_params()
        

    def fill_params(self):
        pars = {'rhox':self._wimp_density,'vE':self._vE,'v0':self._v0,'vesc':self._vesc}
        self._velocity.set_params(pars)

    def set_params(self,pars):
        if 'vE' in pars.keys():
            self._vE = pars['vE']
        if 'v0' in pars.keys():
            self._v0 = pars['v0']
        if 'vesc' in pars.keys():
            self._vesc = pars['vesc']
        if 'rhox' in pars.keys():
            self._wimp_density = pars['rhox']
        self._velocity.set_params(pars)

    def vE(self):
        return self._vE
 
    def v0(self):
        return self._v0

    def vesc(self)
        return self._vesc

    def velocity(self)
        return self._velocity
