""" SolidSphereFormFactor.py

The form factor of a solid sphere charge distribution.

"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from . import FormFactor
from .. import units

class SolidSphereFormFactor(FormFactor):
    """
    This returns the form factor of a solid sphere.

    The solid sphere form factor has the form
    |F(Q^2)|^2 = (j1(Qrn) / Qrn )^2 

    where rn takes the form rn = a A^(1/3) + b
    and represents the atomic radius.

    Attributes:
        A: The mass number
        a: Used in radius calculation
        b: Used in radius calculation
     
    """

    def __init__(self):
        """Initialize to hydrogen with a radius of
           1 fm.
        """

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
                SolidSphereFFA: a 
                SolidSphereFFB: b
        """
        if 'AtomicNumber' in pars.keys():
            self.A = pars['AtomicNumber']
        if 'SolidSphereFFA' in pars.keys():
            self.a = pars['SolidSphereFFA']
        if 'SolidSphereFFB' in pars.keys():
            self.b = pars['SolidSphereFFB']
        self.calculate_params()

    def calculate_params(self):
        """ Calculate the radius given current
            parameters
        """
        self.rn = (self.a * self.A**(1./3) + self.b)

    def initialize(self):
        """ Same as calculate_params """
        self.calculate_params()

    def ff2(self,Q2):
        """ Calculate |F(Q^2)|^2.
 
            Args:
                Q2: Recoil squared momentum transfer
        """

        qrn = np.sqrt(Q2)*self.rn / units.hbarc
         
        j1 = np.sin(qrn) / (qrn*qrn) - np.cos(qrn) / qrn

        return 9 * j1*j1 / (qrn*qrn)
