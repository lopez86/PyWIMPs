""" EngelFormFactor.py

A parameterized spin-dependent form factor mentioned by
Lewin and Smith and attributed to Engel.

"""
__author__ = 'Jeremy P. Lopez'
__date__ = 'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'
from . import FormFactor
from .. import units


class EngelFormFactor(FormFactor):
    """
    This returns a spin-dependent form factor.

    The definitions of the parameters are found
    in Lewin and Smith. This is attributed to Engel.

    This form factor is a sinc function with the
    first zero suppressed.

    Attributes:
        A: The mass number
        a: Used in radius calculation
        b: Used in radius calculation
     
    """


    def __init__(self):
        self.A = 1
        self.a = 1.0 * units.fm
        self.b = 0
        self.calculate_params()

    def set_params(self,pars):
        """ Set parameters from a dictionary.

            Args:
                pars: {string}

            Parameters:
                AtomicNumber: A
                EngelFFRnA: a
                EngelFFRnB: b
        """
        if 'AtomicNumber' in pars.keys():
            self.A = pars['AtomicNumber']
        if 'EngelFFRnA' in pars.keys():
            self.a = pars['EngelFFRnA']
        if 'EngelFFRnB' in pars.keys():
            self.b = pars['EngelFFRnB']
        self.calculate_params()

    def initialize(self):
        """ Same as calculate_params()."""
        self.calculate_params()

    def calculate_params(self):
        """ Calculate some useful parameters."""
        self.rn = (self.a * self.A ** (1./3) + self.b) 

    def ff2(self,Q2):
        """ Calculate |F(Q^2)|^2.
 
            Args:
                Q2: Recoil squared momentum transfer
        """

        qrn = np.sqrt(Q2) * rn / units.hbarc
        if 2.55 < qrn < 4.5:
          return 0.047
        j0 = np.sin(qrn) / qrn
        return j0 * j0
