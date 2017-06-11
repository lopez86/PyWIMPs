from . import FormFactor
from .. import units

class EngelFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.calculate_params()
        self.a = 1.0 * units.fm
        self.b = 0

    def set_params(self,pars):
        if 'AtomicNumber' in pars.keys():
            self.A = pars['AtomicNumber']
        if 'EngelFFRnA' in pars.keys():
            self.a = pars['EngelFFRnA']
        if 'EngelFFRnB' in pars.keys():
            self.b = pars['EngelFFRnB']
        self.calculate_params()

    def calculate_params(self):

        self.rn = (self.a * self.A ** (1./3) + self.b) 

    def ff2(self,Q2):
        qrn = np.sqrt(Q2) * rn / units.hbarc
        if 2.55 < qrn < 4.5:
          return 0.047
        j0 = np.sin(qrn) / qrn
        return j0 * j0
