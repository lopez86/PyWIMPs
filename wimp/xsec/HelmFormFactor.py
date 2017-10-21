""" HelmFormFactor.py

The form factor proposed by Helm discussed in the 
Lewin & Smith dark matter review.

"""
__author__ = 'Jeremy P. Lopez'
__date__ = 'June, 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from . import FormFactor
from .. import units
import numpy as np

class HelmFormFactor(FormFactor):
    """
    This returns the Helm form factor.

    The definitions of the parameters are found
    in Lewin and Smith.

    The Helm form factor has the form
    |F(Q^2)|^2 = 9 |j1(Qrn) / Qrn |^2 exp(-(Qs)^2)

    Attributes:
        A: The mass number
        ca: Used in radius calculation
        cb: Used in radius calculation
        a: Used in radius calculation
        s: Used in radius calculation
     
    """

    def __init__(self):
        """ Initialize to hydrogen and the default
            numbers proposed by Lewin and Smith.
        """
        self.A = 1
        self.ca = 1.23 * units.fm
        self.cb = 0.6 * units.fm
        self.a = 0.52 * units.fm
        self.s = 0.9 * units.fm
        self.calculate_params()
        self.random = np.random

    def set_params(self,pars):
        """ Set parameters from a dictionary.

            Args:
                pars: {string}

            Parameters:
                AtomicNumber: A
                HelmCa: ca
                HelmCb: cb
                HelmA: a
                HelmS: s
 
        """
        if 'AtomicNumber' in pars:
            self.A = pars['AtomicNumber']
        if 'HelmCa' in pars:
            self.ca = pars['HelmCa']
        if 'HelmCb' in pars:
            self.cb = pars['HelmCb']
        if 'HelmA' in pars:
            self.a = pars['HelmA']
        if 'HelmS' in pars:
            self.s = pars['HelmS']

        self.calculate_params()
    
    def initialize(self):
        """ Same as calculate_params()"""
        self.calculate_params()

    def calculate_params(self):
        """ Calculate some useful parameters. """
        self.c = (self.ca * self.A ** (1./3)  - self.cb) 
        self.rn = np.sqrt(self.c*self.c + 7./3 * \
                  np.pi * np.pi * self.a * self.a - \
                  5 * self.s * self.s)
        self.rms = np.sqrt(0.6 * self.rn*self.rn + 3 * self.s * self.s )


    def ff2(self,Q2):
        """ Calculate |F(Q^2)|^2.
 
            Args:
                Q2: Recoil squared momentum transfer
        """
        qrn = np.sqrt(Q2)*self.rn / units.hbarc
         
        j1 = np.sin(qrn) / (qrn*qrn) - np.cos(qrn) / qrn

        return 9 * j1*j1 / (qrn*qrn) * np.exp(-Q2 * self.s*self.s / (units.hbarc*units.hbarc))
