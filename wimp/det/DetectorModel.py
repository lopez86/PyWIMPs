from Efficiency import Efficiency
from Response import Response
from ..mc.sample import Sample
import numpy as np

class DetectorModel:
    def __init__(self):
        self.efficiency = Efficiency()
        self.response = Response()

    def weighted_throw(self,sample):
        sample.weight = self.Efficiency.get_eff(sample) * sample.weight
        return self.response.throw(sample)

    def unweighted_throw(self,sample):
        eff = self.Efficiency.get_eff(sample)
        rnd = np.random.rand()
        if eff < rnd:
            sample.weight = 0
        return self.response.throw(sample)
