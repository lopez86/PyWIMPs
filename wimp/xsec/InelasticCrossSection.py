""" InelasticCrossSection.py

Basic inelastic cross section model: Isotropic in center of mass frame
with a small mass splitting between two dark matter states
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

from . import CrossSection
from .. import units
import numpy as np

class InelasticCS(CrossSection):

    def __init__(self):
        super(InelasticCS,self).__init__()
        self.deltaM = 10*units.keV

    def set_params(self,pars):
        super(InelasticCS,self).set_params(pars)

        if 'MassSplitting' in pars.keys():
            self.deltaM = pars['MassSplitting']

    @property
    def MxPrime(self):
        return self.Mx + self.deltaM

    @property
    def PthresholdCM(self):
        return np.sqrt(2*self.MuPrime * self.deltaM)

    @property
    def EthresholdCM(self):
        return self.MuPrime * self.deltaM / self.Mx

    @property
    def VthresholdCM(self)
        return self.PthresholdCM / self.Mx

    @property
    def Pthreshold(self):
        return self.Mx *np.sqrt(2*self.deltaM/self.Mu)

    @property
    def Ethreshold(self):
        return self.Mx * self.deltaM /self.Mu

    @property
    def Vthreshold(self)
        return units.speed_of_light 
               * np.sqrt(2*self.deltaM/self.Mu)

    @property
    def Mu(self):
        return self.Mx*self.Mt/(self.Mx+self.Mt)

    @property
    def MuPrime(self): 
        return self.MxPrime*self.Mt/(self.MxPrime+self.Mt);


    def MinEr(self,Ex):       
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime-u*u+2*u*vprime)


    def MaxEr(self,Ex):       
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime+u*u+2*u*vprime)


    def dSdErLab(self,Ex,Er):
        if Ex <= self.Ethreshold: return 0
        if Er < minE or Er > maxE: return 0
        minE = self.MinEr(Ex)
        maxE = self.MaxEr(Ex)        
        return self.totalxs / (maxE - minE)

    def PfCM(self,Ex):
        return np.sqrt(2*self.Mx*Ex*self.Mu*self.MuPrime
                       -2*self.MuPrime*self.deltaM)   

    def PiCM(self,Ex):
        return self.Mu * np.sqrt(2*self.Mx*Ex)

    def ErLab(self,Ex,cosTh):
        print("InelasticCS::ErLab(Ex,cosTh_lab) not yet implemented\n")
        return -1

    def cosThetaCMFromCosTheta(self,Ex,cosTh)
        print("InelasticCS::cosThetaCMFromCosTheta not yet implemented\n")
        return -2

    def cosThetaLab(self,Ex,Er):
        cosThCM = self.cosThetaCMFromEr(Ex,Er)
        return self.cosThetaFromCosThetaCM(Ex,cosThCM)

    def cosThetaCMFromEr(self,Ex,Er):
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return (2*Er/self.Mt - vprime*vprime-u*u)/(2*vprime*u)

    def ErFromCosThetaCM(self,Ex,cosTh):
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime+u*u+2*u*vprime*cosTh)

    def cosThetaFromCosThetaCM(self,Ex,cosTh):
        pf = self.PfCM(Ex)
        pi = self.PiCM(Ex)
        return (pf*cosTh+pi) / np.sqrt(pf*pf+pi*pi+2*pf*pi*cosTh)
