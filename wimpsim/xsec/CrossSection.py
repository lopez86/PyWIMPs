import numpy as np
from .. import units

class CrossSection:
    """
    This is the base class for a cross section from a fundamental
    interaction.
    The base class assumes an interaction cross section that is 
    isotropic in the lab frame. For WIMP scattering, this is 
    the typical cross section at 0 momentum transfer used in
    calculations. For spin-1/2 WIMPs scattering on quarks or
    nucleons, all but the pseudoscalar terms have an isotropic
    contribution in the nonrelativistic limit.

    """

    def __init__(self):
      self.totalxs = 1 #Unit normalization is default
      self.Mt = 100*units.GeV
      self.Mx = 100*units.GeV

    def set_params(self,pars):
      if 'XS' in pars.keys():
        self.totalxs = pars['XS']
      if 'Mt' in pars.keys():
        self.Mt = pars['Mt']
      if 'Mx' in pars.keys():
        self.Mx = pars['Mx']

    def dSdOmegaCM(self):
      return self.totalxs / (4*np.pi)

    def dSdcosThCM(self):
      return 0.5*self.totalxs
 

    def MaxEr(self,Ex):
      return 4 * Ex * self.Mt * self.Mx / (self.Mt+self.Mx)**2

    def dSdcosThLab(self,cosTh):
      return 4 * self.totalxs * cosTh

    def dSdErLab(self,Ex,Er):
      maxE = self.MaxEr(Ex)
      if Er >= maxE or Er <0:
        return 0
      return self.totalxs / maxE
 
    def cosThetaLab(self,Ex,Er):
      return math.sqrt(0.5 * Er / self.MaxEr(Ex) )

    def ErLab(self,Ex,cosTh):
      return 2 * self.MaxEr(Ex) * cosTh * cosTh

    def cosThetaCMFromEr(self,Ex,Er):
      return 2 * Ex / self.MaxEr(Ex) - 1

    def cosThetaCMFromCosTheta(self,cosTh):
      return 4 * cosTh*cosTh - 1

    def ErFromCosThetaCM(self,Ex,cosTh):
      return 0.5*self.MaxEr(Ex) * (1+cosTh)   

    def cosThetaFromCosThetaCM(self,cosTh):
      return 0.5 * math.sqrt( 1 + cosTh)
