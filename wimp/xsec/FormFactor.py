""" FormFactor.py

    Form factor calculations.
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np 


class FormFactor:
    """
    Base class for form factor calculations. 

    This is the simplest case, where the form factor is always 1.

    """

    def __init__(self):
        """ Initialize. Doesn't do much here"""
        self._rand = np.random
 
    def initialize(self):
        """ Do any initial calculations of parameters."""
        pass

    @property
    def random(self):
        """ Random number generator."""
        return self._rand

    @random.setter
    def set_random(self,r):
        """ Set the random number generator.

            Args:
                r: (Numpy RandomState)
        """
        self._rand = r

    def set_params(self,pars):
        """ Set parameters with a dictionary.

            Args:
                pars: {string}
        """

        pass

    def ff2(self,Q2):
        """
        Calculates the form factor at the given value of Q^2.

        Q^2 is the squared momentum transfer,
        Q^2 = -q^2 = -(pf - pi)^2 = -2M^2 + 2pf*pi = -2M^2 + 2M(M+E)
        so, Q^2 = 2ME where M is the weight of the nucleus and E is 
        the recoil energy


        Args: 
        Q^2: Squared momentum transfer.
        
        Returns:
        The form factor at Q^2
        """
        return 1
