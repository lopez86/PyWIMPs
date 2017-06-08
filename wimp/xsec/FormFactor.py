class FormFactor:
    """
    Base class for form factor calculations. 

    This is the simplest case, where the form factor is always 1.

    """

    def __init__(self):
        """Doesn't do anything for now."""
        pass

    def set_params(self,pars):
        pass

    def ff2(self,Q2):
        """
        Calculates the form factor at the given value of Q^2.
        Q^2 is the squared momentum transfer,
        Q^2 = -q^2 = -(pf - pi)^2 = -2M^2 + 2pf*pi = -2M^2 + 2M(M+E)
        so, Q^2 = 2ME where M is the weight of the nucleus and E is 
        the recoil energy


        Arguments: 
        Q^2: Squared momentum transfer.
        
        Returns:
        The form factor at Q^2
        """
        return 1
