""" CrossSection.py

    Standard cross section model: Isotropic in the center of mass
    frame at 0 momentum transfer 
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"


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

    Note: In this class, the center of mass angle is the angle
    between the incoming WIMP direction and the outgoing nuclear
    recoil direction.

    """

    def __init__(self):
        """ Initialize to some default values:
            100 GeV masses, cross section of 1 cm^2
        """
        self._rand = np.random
        self.totalxs = 1 #Unit normalization is default
        self.Mt = 100*units.GeV
        self.Mx = 100*units.GeV

    @property
    def random(self):
        """ Random number generator. Not used in base class."""
        return self._rand 

    @random.setter
    def set_random(self,r):
        """ Set the random number generator.
    
            Args:
                r: (Numpy RandomState)
        """
        self._rand = r

    def set_params(self,pars):
        """ Set parameters using a dictionary.

            Args:
                pars: {string}

            Parameters:
                XS: Total cross section
                Mt: Target nucleus mass
                Mx: WIMP mass
        """
        if 'XS' in pars.keys():
          self.totalxs = pars['XS']
        if 'Mt' in pars.keys():
          self.Mt = pars['Mt']
        if 'Mx' in pars.keys():
          self.Mx = pars['Mx']

    def dSdOmegaCM(self):
        """ d(sigma)/d(Omega_CM)

            Double differential cross section in the center
            of mass frame.
        """
        return self.totalxs / (4*np.pi)

    def dSdcosThCM(self):
        """ d(sigma)/d(sin(theta_CM))

            Differential cross section in the center of mass
            frame. 
        """
        return 0.5*self.totalxs
 

    def MaxEr(self,Ex):
        """ Maximum nuclear recoil energy.

            Args:
                Ex: WIMP kinetic energy
        """
        return 4 * Ex * self.Mt * self.Mx / (self.Mt+self.Mx)**2

    def dSdcosThLab(self,cosTh):
        """ Diff. cross section w.r.t. lab frame angle.

            Args:
                cosTh: Center of mass frame cos(theta)
        """
        return 4 * self.totalxs * cosTh

    def dSdErLab(self,Ex,Er):
        """ Diff. cross section w.r.t. nuclear recoil
            kinetic energy.

            Args:
                Ex: WIMP kinetic energy
                Er: Nuclear recoil kinetic energy
        """
        maxE = self.MaxEr(Ex)
        if Er >= maxE or Er <0:
          return 0
        return self.totalxs / maxE
 
    def cosThetaLab(self,Ex,Er):
        """ Lab frame angle.

            Args:
                Ex: WIMP kinetic energy 
                Er: Nuclear recoil kinetic energy
        """
        return np.sqrt(0.5 * Er / self.MaxEr(Ex) )

    def ErLab(self,Ex,cosTh):
        """ Nuclear recoil kinetic energy

            Args:
                Ex: WIMP kinetic energy
                cosTh: Lab frame angle - cos(theta)
        """
        return 2 * self.MaxEr(Ex) * cosTh * cosTh

    def cosThetaCMFromEr(self,Ex,Er):
        """ Center of mass angle

            See note in class documentation on the
            definition of the center of mass angle.

            Args:
                Ex: WIMP kinetic energy
                Er: Nuclear recoil kinetic energy

        """
        return 2 * Ex / self.MaxEr(Ex) - 1

    def cosThetaCMFromCosTheta(self,cosTh):
        """ Center of mass angle

            See note in class documentation on the
            definition of the center of mass angle.

            Args:
                cosTh: Lab frame angle [cos(theta)]
        """
        return 4 * cosTh*cosTh - 1

    def ErFromCosThetaCM(self,Ex,cosTh):
        """ Nuclear recoil kinetic energy

            Args:
                Ex: WIMP kinetic energy
                cosTh: Center of mass angle [cos(theta)]

        """
        return 0.5*self.MaxEr(Ex) * (1+cosTh)   

    def cosThetaFromCosThetaCM(self,cosTh):
        """ Lab frame recoil angle

            Args: 
                cosTh: Center of mass angle [cos(theta)] 
        """
        return 0.5 * np.sqrt( 1 + cosTh)
