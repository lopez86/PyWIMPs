import wimpsim.xsec.FormFactor as FormFactor
import numpy as np

class ExpoFormFactor(FormFactor):
    """
    This is a simple form factor that follows an exponential form.

    |F(Q^2)|^2 = exp(- 2Q^2 / A^2)
    where A is the scale.

    """
    def __init__(self):
        self.scale = 1.0 * units.GeV

    def set_params(self,pars):
        if 'scale' in pars.keys():
            self.scale = pars['scale']

    def ff2(self,Q2):
        return np.exp(-2 * Q2 / (scale * scale))    
