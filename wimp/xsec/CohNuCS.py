""" CohNuCS.py

    Coherent neutrino elastic scattering cross section.

    Only the differential cross section with respect to recoil energy
    is implemented. Beyond this, there are functions to get the lab
    frame recoil angle and the maximum recoil energy.

    Without a form factor, the differential cross section is 
    quadratic in recoil energy, so drawing from this distribution is 
    quite easy.

    Based on:
    Journal of Physics: Conference Series 606 (2015) 012010
"""

__author__    = "Jeremy P. Lopez"
__date__      = "October 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"


import numpy as np
from .. import units
from . import CrossSection

class CohNuCS(CrossSection):
    """
    This is the cross section for neutrinos scattering coherently
    off of a nucleus.
    """

    def __init__(self):
        """ Initialize to some default values:
            100 GeV masses, cross section of 1 cm^2
        """

        super(CohNuCS,self).__init__()
        self.Qw = 1
        self.GF = 1.1663787e-5 / (units.GeV*units.Gev)
                  * units.hbar*units.hbar

    def initialize(self):
        """ Do any initial calculations of parameters."""
        pass

    @property
    def random(self):
        """ Random number generator. Not used in base class."""
        return self._rand 

    @random.setter
    def random(self,r):
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
                XS: Total cross section (has no effect)
                Mt: Target nucleus mass
                Mx: WIMP mass (has no effect)
                Qw: Nuclear weak charge
                
        """
        super(CohNuCS,self).set_params(pars)

        if 'Qw' in pars:
            self.Qw = pars['Qw']
        
    def dSdOmegaCM(self):
        """ d(sigma)/d(Omega_CM). Not implemented.

            Double differential cross section in the center
            of mass frame.
        """
        print('dSdOmegaCM not implemeneted for CEvNS')
        return -1.0

    def dSdcosThCM(self):
        """ d(sigma)/d(sin(theta_CM)). Not implemented.

            Differential cross section in the center of mass
            frame. 
        """
        print('dSdcosThCM not implemented for CEvNS')
        return -1.0
 

    def MaxEr(self,Ev):
        """ Maximum nuclear recoil energy.

            Args:
                Ev: Neutrino kinetic energy
        """
        return 2*Ev*Ev/self.Mt

    def dSdcosThLab(self,cosTh):
        """ Diff. cross section w.r.t. lab frame angle. Not implemented.

            Args:
                cosTh: Center of mass frame cos(theta)
        """
        print('dSdcosThLab not implemented for CEvNS')
        return -1.0

    def dSdErLab(self,Ev,Er):
        """ Diff. cross section w.r.t. nuclear recoil
            kinetic energy.

            Args:
                Ev: Neutrino kinetic energy
                Er: Nuclear recoil kinetic energy
        """
        maxE = self.MaxEr(Ev)
        if Er >= maxE or Er <0:
          return 0
        # Leading constant
        const = 0.125/np.pi * self.GF*self.GF*self.Mt \
                * self.Qw*self.Qw
        # Kinematic part
        kin = 2 - 2*Er/Ev + Er*Er/(Ev*Ev) -self.Mt * Er/(Ev*Ev)

        # Not included: form factor

        return const * kin
 
    def cosThetaLab(self,Ev,Er):
        """ Lab frame angle.

            Args:
                Ev: Neutrino kinetic energy 
                Er: Nuclear recoil kinetic energy
        """
        pr = np.sqrt(2*Er*self.Mt)
        return (2*(self.Mt+Ev)*Er - Er*Er) / (2*Ev*pr)

    def ErLab(self,Ev,cosTh):
        """ Nuclear recoil kinetic energy. Not Implemented.

            Args:
                Ex: Neutrino kinetic energy
                cosTh: Lab frame angle - cos(theta)
        """
        print('ErLab not implemented for CEvNS')
        return -1

    def cosThetaCMFromEr(self,Ev,Er):
        """ Center of mass angle. Not implemented.

            See note in class documentation on the
            definition of the center of mass angle.

            Args:
                Ex: Neutrino kinetic energy
                Er: Nuclear recoil kinetic energy

        """
        
        print('cosThetaCMFromEr not implemented for CEvNS')
        return -2.0

    def cosThetaCMFromCosTheta(self,cosTh):
        """ Center of mass angle. Not implemented.

            See note in class documentation on the
            definition of the center of mass angle.

            Args:
                cosTh: Lab frame angle [cos(theta)]
        """
        
        print('cosThetaCMFromCosTheta not implemented for CEvNS')
        return -2.0

    def ErFromCosThetaCM(self,Ev,cosTh):
        """ Nuclear recoil kinetic energy. Not implemented.

            Args:
                Ex: Neutrino kinetic energy
                cosTh: Center of mass angle [cos(theta)]

        """

        print('ErFromCosThetaCM not implemented for CEvNS')
        return -1.0

    def cosThetaFromCosThetaCM(self,cosTh):
        """ Lab frame recoil angle. Not implemented.

            Args: 
                cosTh: Center of mass angle [cos(theta)] 
        """
        
        print('cosThetaFromCosThetaCM not implemented for CEvNS')
        return -2.0
