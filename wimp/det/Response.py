import numpy as np
from ..mc.sample import Sample

class Response:

    def __init__(self):
        pass

    def throw(self,s):
        return s


class GaussianResponse:
    def __init__(self):
        self.sigma = 0.1

    def throw(self,s):
        s.Er = np.random.normal(s.Er,self.sigma * s.Er)
        return s
