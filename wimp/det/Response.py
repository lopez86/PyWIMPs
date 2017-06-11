import numpy as np
from ..mc.sample import Sample

class Response:

    def __init__(self):
        pass

    def set_params(self,pars):
        pass

    def throw(self,s):
        return s


class GaussianResponse:
    def __init__(self):
        self.sigma = 0.1

    def set_params(self,pars):
        if 'RespSigma' in pars.keys():
            self.sigma = pars['RespSigma']

    def throw(self,s):
        s.Er = np.random.normal(s.Er,self.sigma * s.Er)
        return s
