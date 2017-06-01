from .. import units
from .CrossSection import CrossSection
from .FormFactor import FormFactor

class InteractionModel:

    def __init__(self):
        self._cross_section = CrossSection()
        self._form_factor = FormFactor()
        self._Mx = 100*GeV
        self._Mt = 100*GeV
        self._Mtot = 100 * kg
        self._total_xs = 1e-40*cm*cm
        self.fill_params()

    def fill_params(self):
       pars = {"Mx":self._Mx,"Mt":self._Mt,"XS":self._total_xs }
       self._cross_section.set_params(pars)
       self._form_factor.set_params(pars)


    def Mx(self): 
        return self._Mx
    def Mt(self): 
        return self._Mt
    def Mtot(self): 
        return self._Mtot
    def total_xs(self): 
        return self._total_xs
    def form_factor(self): 
        return self._form_factor
    def cross_section(self):
        return self._cross_section
    
    def set_Mx(self,m):
        self._Mx = m
        par = {"Mx":m,}
        self._cross_section().set_params(par)
        self._form_factor().set_params(par)

    def set_Mt(self,m):
        self._Mt = m
        par = {"Mt":m,}
        self._cross_section().set_params(par)
        self._form_factor().set_params(par)


    def set_total_xs(self,xs):
        self._total_xs = xs
        par = {"XS":xs,}
        self._cross_section().set_params(par)
        self._form_factor().set_params(par)

    def set_form_factor(self,ff):
        self._form_factor = ff
        self.fill_params()

    def set_cross_section(self,xs):
        self._cross_section = xs
        self.fill_params()
