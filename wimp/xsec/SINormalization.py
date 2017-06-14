""" SINormalization.py

    Normalize a nuclear cross section to a nucleon 
    cross section.
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

from .. import units

class SINormalization:
    """ Class to normalize a spin-independent cross section 
        to a spin-independent nucleon cross section.

        Attributes:
           gp: Proton coupling
           gn: Neutron coupling
    """
    def __init__(self):
        """ Initialize. 
            Default is equal proton/neutron couplings.
        """
        self.gp = 1
        self.gn = 1

    def normalize(self,nucleus,Mx,to_proton=True):
        """ Get the normalization constant.

            Args:
                nucleus: A Nucleus object
                to_proton: If true, normalize to proton
                           else normalize to neutron
        """
        m = nucleus.mass
        mu = m * Mx / (Mx + m)
        muN = Mx * units.proton_mass / (units.proton_mass + Mx)
        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)

        g = nucleus.Z * self.gp + nucleus.N * self.gn
        gN = self.gp

        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)
            gN = self.gn
    
        return muN*muN * gN * gN / ( mu*mu*g*g)
        
