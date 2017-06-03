from . import FormFactor
from .. import units

class EngelFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.calculate_params()
        self.a = 1.0 / units.hbarc
        self.b = 0

    def set_params(self,pars):
        if 'A' in pars.keys():
            self.A = pars['A']
        if 'a' in pars.keys():
            self.a = pars['a']
        if 'b' in pars.keys():
            self.b = pars['b']
        self.calculate_params()

    def calculate_params(self):

        self.rn = self.a * self.A ** (1./3) + self.b

    def ff2(self,Q2):
        qrn = np.sqrt(Q2) * rn
        if 2.55 < qrn < 4.5:
          return 0.047
        j0 = np.sin(qrn) / qrn
        return j0 * j0
