import numpy as np
from ..mc.sample import Sample

class Efficiency:
    def __init__(self):
        self.eff = 1

    def get_eff(self,s):
        return self.eff

class LogisticEfficiency(Efficiency)
    def __init__(self):
        self.eff_max = 1
        self.x0 = 5
        self.xscale=4

    def get_eff(self,s):
        return 1./(1 +np.exp(-(s.Er - self.x0)/self.xscale )
    
