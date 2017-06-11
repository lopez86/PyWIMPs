from . import FormFactor
from .. import units

class ThinShellFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.a = 1.0 * units.fm
        self.b = 0
        self.calculate_params()


    def set_params(self,pars):
        if 'AtomicNumber' in pars.keys():
            self.A = pars['AtomicNumber']
        if 'ThinShellFFA' in pars.keys():
            self.a = pars['ThinShellFFA']
        if 'ThinShellFFB' in pars.keys():
            self.b = pars['ThinShellFFB']
        self.calculate_params()

    def calculate_params(self):

        self.rn = self.a * self.A**(1./3) + self.b



    def ff2(self,Q2):
        qrn = np.sqrt(Q2)*self.rn / units.hbarc
         
        f = np.sin(qrn) / qrn

        return f*f
