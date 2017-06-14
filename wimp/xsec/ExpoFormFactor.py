""" ExpoFormFactor.py

A form factor that is exponential in Q^2.
"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from .FormFactor import FormFactor
import numpy as np


class ExpoFormFactor(FormFactor):
    """
    This is a simple form factor that follows an exponential form.

    |F(Q^2)|^2 = exp(- 2Q^2 / A^2)
    where A is the scale.

    Attributes:
        scale: The scale in energy units.

    """
    def __init__(self):
        """ Initialize to a scale of 1 GeV"""
        self.scale = 1.0 * units.GeV

    def set_params(self,pars):
        """ Set parameters from a dictionary.

            Args:
                pars: {string}

            Parameters:
                ExpoFFScale
        """

        if 'ExpoFFScale' in pars.keys():
            self.scale = pars['ExpoFFScale']

    def ff2(self,Q2):
        """ Calculate |F(Q^2)|^2.
 
            Args:
                Q2: Recoil squared momentum transfer
        """

        return np.exp(-2 * Q2 / (scale * scale))    
