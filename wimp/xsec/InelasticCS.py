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
    """ Cross section for a basic inelastic dark matter model.

    This model assumes that there is some small mass splitting 
    between dark matter WIMPs and a second type of WIMP. For 
    conservation of energy, there will be a non-zero velocity 
    threshold for the interaction to happen.

    """
    def __init__(self):
        """
        Initialize to some default values:

        100 GeV masses, 10 keV splitting, 1 cm^2 cross section
        """
        super(InelasticCS,self).__init__()
        self.deltaM = 10*units.keV

    def set_params(self,pars):
        """ Set parameters from a dictionary.
 
        Args:
            pars: {string}

        Parameters:
            MassSplitting: Mass difference between the two states
            Other CrossSection parameters
        """
        super(InelasticCS,self).set_params(pars)

        if 'MassSplitting' in pars:
            self.deltaM = pars['MassSplitting']

    @property
    def MxPrime(self):
        """ Final state WIMP mass"""
        return self.Mx + self.deltaM

    @property
    def PthresholdCM(self):
        """ Momentum threshold in CM frame"""
        return np.sqrt(2*self.MuPrime * self.deltaM)

    @property
    def EthresholdCM(self):
        """ WIMP energy threshold in CM frame"""
        return self.MuPrime * self.deltaM / self.Mx

    @property
    def VthresholdCM(self):
        """ WIMP velocity threshold in CM frame"""
        return self.PthresholdCM / self.Mx

    @property
    def Pthreshold(self):
        """ WIMP momentum threshold in lab frame"""
        return self.Mx *np.sqrt(2*self.deltaM/self.Mu)

    @property
    def Ethreshold(self):
        """ WIMP kinetic energy threshold in lab frame"""
        return self.Mx * self.deltaM /self.Mu

    @property
    def Vthreshold(self):
        """ WIMP velocity threshold in lab frame"""
        return units.speed_of_light \
               * np.sqrt(2*self.deltaM/self.Mu)

    @property
    def Mu(self):
        """ Initial state reduced mass"""
        return self.Mx*self.Mt/(self.Mx+self.Mt)

    @property
    def MuPrime(self): 
        """ Final state reduced mass"""
        return self.MxPrime*self.Mt/(self.MxPrime+self.Mt);


    def MinEr(self,Ex):       
        """ Minimum recoil energy given a WIMP energy

        Args:
            Ex: WIMP kinetic energy
        """
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime-u*u+2*u*vprime)


    def MaxEr(self,Ex):       
        """ Maximum recoil energy given a WIMP energy
 
        Args:
            Ex: WIMP kinetic energy
        """
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime+u*u+2*u*vprime)


    def dSdErLab(self,Ex,Er):
        """ Single differential cross section with respect to
            recoil energy

        Args:
            Ex: WIMP kinetic energy
            Er: WIMP recoil energy
        """
        if Ex <= self.Ethreshold: return 0
        if Er < minE or Er > maxE: return 0
        minE = self.MinEr(Ex)
        maxE = self.MaxEr(Ex)        
        return self.totalxs / (maxE - minE)

    def PfCM(self,Ex):
        """ Center of mass momentum in final state
      
        Args:
            Ex: WIMP kinetic energy
        """
        return np.sqrt(2*self.Mx*Ex*self.Mu*self.MuPrime
                       -2*self.MuPrime*self.deltaM)   

    def PiCM(self,Ex):
        """ Center of mass momentum in initial state
      
        Args:
            Ex: WIMP kinetic energy
        """
        return self.Mu * np.sqrt(2*self.Mx*Ex)

    def ErLab(self,Ex,cosTh):
        """ Not implemented at this point. Recoil energy as 
            function of recoil angle

        Args:
            Ex: WIMP kinetic energy
            cosTh: Recoil cos(theta) in lab frame
        """
        print("InelasticCS::ErLab(Ex,cosTh_lab) not yet "
              "implemented\n")
        return -1

    def cosThetaCMFromCosTheta(self,Ex,cosTh):
        """ Not implemented at this point. CM recoil angle as
            function of recoil angle

        Args:
            Ex: WIMP kinetic energy
            cosTh: Recoil cos(theta) in lab frame
        """
        print("InelasticCS::cosThetaCMFromCosTheta "
              "not yet implemented\n")
        return -2

    def cosThetaLab(self,Ex,Er):
        """ Lab frame recoil angle as function of recoil energy

        Args:
            Ex: WIMP kinetic energy
            Er: Recoil kinetic energy
        """
        cosThCM = self.cosThetaCMFromEr(Ex,Er)
        return self.cosThetaFromCosThetaCM(Ex,cosThCM)

    def cosThetaCMFromEr(self,Ex,Er):
        """ CM recoil angle as function of lab frame recoil energy
 
        Args:
            Ex: WIMP kinetic energy
            Er: Recoil kinetic energy
        """
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return (2*Er/self.Mt - vprime*vprime-u*u)/(2*vprime*u)

    def ErFromCosThetaCM(self,Ex,cosTh):
        """ Lab frame recoil energy as function of CM recoil angle

        Args:
            Ex: WIMP kinetic energy
            cosTh: Recoil cos(theta) in CM frame
        """
        vprime = self.PfCM(Ex) / self.Mt
        u = self.PiCM(Ex)/self.Mt
        return 0.5*self.Mt*(vprime*vprime+u*u+2*u*vprime*cosTh)

    def cosThetaFromCosThetaCM(self,Ex,cosTh):
        """ Lab frame recoil angle as function of CM recoil angle

        Args:
            Ex: WIMP kinetic energy
            cosTh: Recoil cos(theta) in CM frame
        """
        pf = self.PfCM(Ex)
        pi = self.PiCM(Ex)
        return (pf*cosTh+pi) / np.sqrt(pf*pf+pi*pi+2*pf*pi*cosTh)
