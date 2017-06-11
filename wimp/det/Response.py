import numpy as np
from ..mc.sample import Sample

class Response:

    def __init__(self):
        self.rand = np.random
 
    def set_random(self,r):
        self.rand = r

    def set_params(self,pars):
        pass

    def throw(self,s):
        return s


class GaussianResponse(Response):
    def __init__(self):
        self.rand = self.rand
        self.sigma = 0.1

    def set_params(self,pars):
        if 'RespSigma' in pars.keys():
            self.sigma = pars['RespSigma']

    def throw(self,s):
        s.Er = self.rand.normal(s.Er,self.sigma * s.Er)
        return s
