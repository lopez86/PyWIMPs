""" DipoleFormFactor.py

Dipole form factor calculations.

"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from .FormFactor import FormFactor
from .. import units


class DipoleFormFactor(FormFactor):
    """
    This returns the dipole form factor given an interaction Q2
    and the dipole scale.

    The dipole form factor has the form
    |F(Q^2)|^2 = (1 + Q^2/A^2)^-2 
    and is associated with an exponential charge distribution.

    Attributes:
        scale: The energy scale of the form factor.
     
    """
    def  __init__(self):
        """Initialize with a scale of 1 GeV """
        
        self.scale = 1.0 * units.GeV

    def set_params(self,pars):
        """ Set parameters from a dictionary.

            Args:
                pars: {string}

            Parameters:
                DipoleFFScale
        """
        if 'DipoleFFScale' in pars.keys:
          self.scale = pars['DipoleFFScale']

    def ff2(self,Q2):
        """ Calculate |F(Q^2)|^2.
 
            Args:
                Q2: Recoil squared momentum transfer
        """
        return (1 + Q2/(scale*scale))**-2
