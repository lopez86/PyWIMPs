import FormFactor

class DipoleFormFactor(FormFactor):
    """
    This returns the dipole form factor given an interactoin Q2
    and the dipole scale.

    The dipole form factor has the form
    |F(Q^2)|^2 = (1 + Q^2/A^2)^-2 
    and is associated with an exponential charge distribution.

    """
    def  __init__(self,scale):
        self.scale = scale

    def ff2(self,Q2):
      return (1 + Q2/(scale*scale))**-2
