from .. import units
from .CrossSection import CrossSection
from .FormFactor import FormFactor

class InteractionModel:

    def __init__(self):
        self.cross_section = CrossSection()
        self.form_factor = FormFactor()
        self._Mx = 100*units.GeV
        self._Mt = 100*units.GeV
        self._Mtot = 100 * units.kg
        self._total_xs = 1e-40*units.cm*units.cm
        self.fill_params()

    def fill_params(self):
        pars = {"Mx":self._Mx,"Mt":self._Mt,"XS":self._total_xs}
        self.cross_section.set_params(pars)
        self.form_factor.set_params(pars)

    def set_params(self,pars):

        if 'Mx' in pars.keys():
            self._Mx = pars['Mx']
        if 'Mt' in pars.keys():
            self._Mt = pars['Mt']
        if 'Mtot' in pars.keys():
            self._Mtot = pars['Mtot']
        if 'XS' in pars.keys():
            self._total_xs = pars['XS']

        self.cross_section.set_params(pars)
        self.form_factor.set_params(pars)
    def Mx(self): 
        return self._Mx
    def Mt(self): 
        return self._Mt
    def Mtot(self): 
        return self._Mtot
    def total_xs(self): 
        return self._total_xs

    
    def set_Mx(self,m):
        self._Mx = m
        par = {"Mx":m,}
        self.cross_section.set_params(par)
        self.form_factor.set_params(par)

    def set_Mt(self,m):
        self._Mt = m
        par = {"Mt":m,}
        self.cross_section.set_params(par)
        self.form_factor.set_params(par)

    def set_Mtot(self,m):
        self._Mtot = m

    def set_total_xs(self,xs):
        self._total_xs = xs
        par = {"XS":xs,}
        self.cross_section.set_params(par)
        self.form_factor.set_params(par)

    def set_form_factor(self,ff):
        self.form_factor = ff
        self.fill_params()

    def set_cross_section(self,xs):
        self.cross_section = xs
        self.fill_params()
