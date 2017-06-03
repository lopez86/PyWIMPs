from . import FormFactor
from .. import units

class HelmFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.calculate_params()


    def set_params(self,pars):
        if 'A' in pars.keys():
            self.A = pars['A']
        self.calculate_params()

    def calculate_params(self):

        self.c = (1.23 * self.A ** (1./3)  - 0.6) / units.hbarc
        self.a = 0.52 / units.hbarc
        self.s = 0.9 / hbarc
        self.rn = np.sqrt(c*c + 7./3 * np.pi * np.pi * a * a - 5 * s * s)
        self.rms = 0.6 * rn*rn + 3 * s * s


    def ff2(self,Q2):
        qrn = np.sqrt(Q2)*self.rn
         
        j1 = np.sin(qrn) / (qrn*qrn) - np.cos(qrn) / qrn

        return 9 * j1*j1 / (qrn*qrn) * np.exp(-Q2 * self.s*self.s)
