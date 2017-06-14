""" SDNormalization.py

    Normalize a nuclear cross section to a nucleon 
    cross section.
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np

class SDNormalization:
    """ Class to normalize a spin-dependent cross section
        to a spin-dependent nucleon cross section.

        Several methods are provided, as there are several
        ways described in the literature to do this.
    """    
    def __init__(self):
        """ Initialize. 

            Default is equal proton and neutron couplings.
        """
        self.gp = 1
        self.gn = 1

    def single_particle(self,nucleus,Mx,to_proton=True):
        """ Normalization for single particle model.

            An odd-even nucleus where we assume that 
            J = L_unpaired + S_unpaired

            Args:
                nucleus: Nucleus object
                Mx:      WIMP mass
                to_proton: Normalize to proton if true
                           otherwise to neutron
        """
        m = nucleus.mass
        mu = m * Mx / (Mx + m)
        muN = Mx * units.proton_mass / (units.proton_mass + Mx)
        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)
 
        # Odd type:
        g2 = self.gp*self.gp
        if nucleus.N%2==1:
          g2 = self.gn*self.gn

        g2 *= ( nucleus.J*(nucleus.J+1) + 0.75 -\
                nucleus.L*(nucleus.L+1))**2 /\
               (4*nucleus.J*(nucleus.J+1) )

        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        
    def odd_group(self,nucleus,Mx,to_proton=True):
        """ Normalization for odd-group model.

            An odd-even nucleus where we take the average
            spin from all of the odd group nuclei after
            subtracting angular momentum contributions from
            both the nucleus and the nucleon.

            Args:
                nucleus: Nucleus object
                Mx:      WIMP mass
                to_proton: Normalize to proton if true
                           otherwise to neutron
        """
        m = nucleus.mass
        mu = m * Mx / (Mx + m)
        muN = Mx * units.proton_mass / (units.proton_mass + Mx)
        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)

        # Odd type:
        g2 = self.gp*self.gp * (nucleus.J+1)/nucleus.J
        g2 *= ((nucleus.magnetic_moment - nucleus.J)/(units.proton_mag_moment - 1))**2
        if n.N()%2==1:
            g2 = self.gn*self.gn * (nucleus.J+1)/nucleus.J
            g2 *= (nucleus.magnetic_moment/neutron_mag_moment)**2
 
        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        
    def spin_content(self,nucleus,Mx,to_proton=True):
        """ Normalization for the precise treatment using
            both proton and neutron spin contributions.

            Args:
                nucleus: Nucleus object
                Mx:      WIMP mass
                to_proton: Normalize to proton if true
                           otherwise to neutron
        """
        m = nucleus.mass
        mu = m * Mx / (Mx + m)
        muN = Mx * units.proton_mass / (units.proton_mass + Mx)
        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)


        g2 = (nucleus.J+1) / nucleus.J
        g2 *= (self.gp * nucleus.proton_spin + \
               self.gn*nucleus.neutron_spin)**2

        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * units.neutron_mass / (units.neutron_mass + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        

    def normalization(self,nucleus,Mx,method='spin_content',to_proton=True):
        """ Normalize according to the desired method.

            Args:
                nucleus: Nucleus object
                Mx:      WIMP mass
                method: 'odd_group', 'single_particle', or
                        [DEFAULT] 'spin_content' 
                to_proton: Normalize to proton if true
                           otherwise to neutron
        """
        if method== 'odd_group':
            return self.odd_group(nucleus,Mx,to_proton)
        elif method == 'single_particle':
            return self.single_particle(nucleus,Mx,to_proton)

        return self.spin_content(nucleus,Mx,to_proton)
