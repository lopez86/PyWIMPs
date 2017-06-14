""" InteractionModel.py

    Container for interaction properties.
"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

from .. import units
from .CrossSection import CrossSection
from .FormFactor import FormFactor
import numpy as np

class InteractionModel:
    """ Container class for interaction properties.
 
        Attributes:
            cross_section: CrossSection
            form_factor:   FormFactor
            Mtot:          detector mass
    """

    def __init__(self):
        """ Initialize to some default values.

            100 GeV masses, 100 kg detector,
            10^-40 cm^2 cross section.
        """
        self._rand = np.random
        self.cross_section = CrossSection()
        self.form_factor = FormFactor()
        self._Mx = 100*units.GeV
        self._Mt = 100*units.GeV
        self._Mtot = 100 * units.kg
        self._total_xs = 1e-40*units.cm*units.cm
  
        self.cross_section.Mx = self._Mx
        self.cross_section.Mt = self._Mt
        self.cross_section.totalxs = self._total_xs

    def initialize(self):
        """ Do any initial calculations of parameters."""
        self.cross_section.initialize()
        self.form_factor.initialize()


    @property
    def random(self):
        """ Random number generator. Not used in base class."""
        return self._rand

    @random.setter
    def set_random(self,r):
        """ Set the random number generator.

            Args:
                r: (Numpy RandomState)
        """
        self._rand = r
        self.cross_section.set_random(r)
        self.form_factor.set_random(r)
 
    def set_params(self,pars):
        """ Set parameters using a dictionary.

            Args:
                pars: {string}

            Parameters:
                Mx: WIMP mass
                Mt: Target nucleus mass
                Mtot: Total detector fiducial mass
                XS: Total cross section
        """
        if 'Mx' in pars.keys():
            self.Mx = pars['Mx']
        if 'Mt' in pars.keys():
            self.Mt = pars['Mt']
        if 'Mtot' in pars.keys():
            self.Mtot = pars['Mtot']
        if 'XS' in pars.keys():
            self.total_xs = pars['XS']

        self.cross_section.set_params(pars)
        self.form_factor.set_params(pars)

    @property
    def Mx(self): 
        """ WIMP mass"""
        return self._Mx
   
    @property
    def Mt(self): 
        """ Target nucleus mass"""
        return self._Mt

    @property
    def Mtot(self): 
        """ Detector fiducial mass"""
        return self._Mtot

    @property
    def total_xs(self): 
        """ Total cross section"""
        return self._total_xs

    @Mx.setter
    def set_Mx(self,m):
        """ Set the WIMP mass

            Args:
                m: the mass
        """ 
        self._Mx = m       
        self.cross_section.Mx = m

    @Mt.setter
    def set_Mt(self,m):
        """ Set the target nucleus mass

            Args:
                m: the mass
        """ 
        self._Mt = m
        self.cross_section.Mt = m

    @Mtot.setter
    def set_Mtot(self,m):
        """ Set the detector fiducial mass

            Args:
                m: the mass
        """ 

        self._Mtot = m

    @total_xs.setter
    def set_total_xs(self,xs):
        """ Set the total cross section

            Args:
                xs: the cross section
        """ 

        self._total_xs = xs
        self.cross_section.totalxs = xs

    def set_form_factor(self,ff):
        """ Set the form factor object

            Args:
                ff: FormFactor
        """ 

        self.form_factor = ff

    def set_cross_section(self,xs):
        """ Set the cross section object

            Args:
                xs: CrossSection
        """ 

        self.cross_section = xs
