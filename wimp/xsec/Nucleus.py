""" Nucleus.py

Information about nuclei relevant for dark matter 
calculations.

"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

class Nucleus:
    """ Class to hold information about nuclei.

        Attributes:
            A: Mass number
            Z: Atomic number
            mas: Atomic mass
            magnetic_moment: nuclear magnetic moment
            spin: nuclear spin (total angular momentum)
            Sp: proton spin contribution to spin
            Sn: neutron spin contribution to spin
            L: angular momentum (typically used in single
               particle spin model)
    """
    
    def __init__(self,pars=None):
    """ Initialize the object.

        Args:
            pars: {string}
    """
        if pars is not None:
            self.set_params(pars)
        else:
            self.A = 0
            self.Z = 0
            self.mass = 0
            self.magnetic_moment = 0
            self.spin = 0
            self.Sp = 0
            self.Sn = 0
            self.L = 0
    

    def set_params(self,pars):
        """ Set parameters based on a dictionary.
 
            Args:
                pars: {string}
 
            Parameters:

                NuclMassNumber:   A
                NuclAtomicNumber: Z
                NuclMass:         mass
                NuclMagMoment:    magnetic_moment
                NuclSpin:         spin
                NuclSp:           Sp
                NuclSn:           Sn
                NuclAngMom:       L
        """
        if 'NuclMassNumber' in pars:
            self.A = pars['NuclMassNumber']
        if 'NuclAtomicNumber' in pars:
            self.Z = pars['NuclAtomicNumber']
        if 'NuclMass' in pars:
            self.mass = pars['NuclMass']
        if 'NuclSpin' in pars:
            self.spin = pars['NuclSpin']
        if 'NuclSp' in pars:
            self.Sp = pars['NuclSp']
        if 'NuclSn' in pars:
            self.Sn = pars['NuclSn']
        if 'NuclAngMom' in pars:
            self.L = pars['NuclAngMom']
        if 'NuclMagMoment' in pars:
            self.magnetic_moment = pars['NuclMagMoment']

    @property
    def N(self):
        """ The number of neutrons."""
        return self.A - self.Z

    @property
    def J(self)
        """ The nuclear spin."""
        return self.spin

