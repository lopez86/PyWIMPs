from . import FormFactor
from .. import units

class SolidSphereFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.a = 1.0 / units.hbarc
        self.b = 0
        self.calculate_params()


    def set_params(self,pars):
        if 'A' in pars.keys():
            self.A = pars['A']
        if 'a' in pars.keys():
            self.a = pars['a']
        if 'b' in pars.keys()
            self.b = pars['b']
        self.calculate_params()

    def calculate_params(self):

        self.rn = self.a * self.A**(1./3) + self.b



    def ff2(self,Q2):
        qrn = np.sqrt(Q2)*self.rn
         
        j1 = np.sin(qrn) / (qrn*qrn) - np.cos(qrn) / qrn

        return 9 * j1*j1 / (qrn*qrn)
