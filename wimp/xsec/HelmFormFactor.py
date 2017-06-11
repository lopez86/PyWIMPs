from . import FormFactor
from .. import units

class HelmFormFactor(FormFactor):

    def __init__(self):
        self.A = 1
        self.ca = 1.23 * units.fm
        self.cb = 0.6 * units.fm
        self.a = 0.52 * units.fm
        self.s = 0.9 * units.fm
        self.calculate_params()


    def set_params(self,pars):
        if 'AtomicNumber' in pars.keys():
            self.A = pars['AtomicNumber']
        if 'HelmCa' in pars.keys():
            self.ca = pars['HelmCa']
        if 'HelmCb' in pars.keys():
            self.cb = pars['HelmCb']
        if 'HelmA' in pars.keys():
            self.a = pars['HelmA']
        if 'HelmS' in pars.keys():
            self.s = pars['HelmS']

        self.calculate_params()
    

    def calculate_params(self):

        self.c = (self.ca * self.A ** (1./3)  - self.cb) 
        self.rn = np.sqrt(self.c*self.c + 7./3 * \
                  np.pi * np.pi * self.a * self.a - \
                  5 * self.s * self.s)
        self.rms = np.sqrt(0.6 * rn*rn + 3 * self.s * self.s )


    def ff2(self,Q2):
        qrn = np.sqrt(Q2)*self.rn / units.hbarc
         
        j1 = np.sin(qrn) / (qrn*qrn) - np.cos(qrn) / qrn

        return 9 * j1*j1 / (qrn*qrn) * np.exp(-Q2 * self.s*self.s / (units.hbarc*units.hbarc))
