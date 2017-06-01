from .FormFactor import FormFactor
from .. import units

class DipoleFormFactor(FormFactor):
    """
    This returns the dipole form factor given an interactoin Q2
    and the dipole scale.

    The dipole form factor has the form
    |F(Q^2)|^2 = (1 + Q^2/A^2)^-2 
    and is associated with an exponential charge distribution.

    """
    def  __init__(self):
        self.scale = 1.0 * units.GeV

    def set_params(self,pars):
        if 'scale' in pars.keys:
          self.scale = pars['scale']

    def ff2(self,Q2):
        return (1 + Q2/(scale*scale))**-2
