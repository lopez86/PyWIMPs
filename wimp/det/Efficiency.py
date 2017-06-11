import numpy as np
from ..mc.sample import Sample

class Efficiency:
    def __init__(self):
        self.rand = np.random
        self.eff = 1

    def set_random(self,r):
        self.rand = r

    def set_params(self,pars):
        if 'EffMax' in pars.keys():
            self.eff = pars['EffMax']

    def get_eff(self,s):
        return self.eff

class LogisticEfficiency(Efficiency):
    def __init__(self):
        self.eff = 1
        self.x0 = 5
        self.xscale=4

    def set_params(self,pars):
        if 'EffMax' in pars.keys():
            self.eff = pars['EffMax']
        if 'EffEr0' in pars.keys():
            self.x0 = pars['EffEr0']
        if 'EffErScale' in pars.keys():
            self.xscale = pars['EffErScale']

    def get_eff(self,s):
        return self.eff/(1 +np.exp(-(s.Er - self.x0)/self.xscale ) )
    
