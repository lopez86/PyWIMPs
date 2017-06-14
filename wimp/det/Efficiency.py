""" Efficiency.py

    Basic classes to represent the detector efficiency.


"""
__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np
from ..mc.sample import Sample

class Efficiency:
    """ Base efficiency class: constant efficiency

        Calculates the efficiency for events in a detector.

        Attributes:
            eff: The efficiency

    """
    def __init__(self):
        """ Initializes to 100% efficiency and default numpy 
            randomization
        """
        self._rand = np.random
        self.eff = 1

    def initialize(self):
        """ Do any initial calculations of parameters."""
        pass


    @property
    def random(self):
        """ The random number generator. """
        return self._rand

    @random.setter
    def set_random(self,r):
        """ Set the random number generator.

            Args:
                r: (Numpy RandomState or equivalent)
        """

        self._rand = r

    def set_params(self,pars):
        """ Set the parameters given a dictionary

            Args:
                pars: ({string})

            Parameters:
                EffMax: the efficiency

        """
        if 'EffMax' in pars.keys():
            self.eff = pars['EffMax']

    def efficiency(self,sample):
        """ Get the efficiency for the given event using a 
            constant value.

            Args:
                sample: (PyWIMPs Sample)

            Returns: the efficiency
        """
        return self.eff

class LogisticEfficiency(Efficiency):
    """ Calculate the efficiency given a logistic model in recoil
        energy.

        Attributes:
            eff: (double) - the maximum efficiency
            x0: (double) - the energy offset
            xscale: (double) - the energy scale
    """
    def __init__(self):
        """ Set some initial values"""
        self.eff = 1
        self.x0 = 5
        self.xscale=4

    def set_params(self,pars):
        """ Set parameters given a dictionary 

            Args:
                pars: {string}

            Parameters:
                EffMax: the maximum efficiency
                EffEr0: the energy offset
                EffErScale: the energy scale
        """
      
        if 'EffMax' in pars.keys():
            self.eff = pars['EffMax']
        if 'EffEr0' in pars.keys():
            self.x0 = pars['EffEr0']
        if 'EffErScale' in pars.keys():
            self.xscale = pars['EffErScale']

    def efficiency(self,sample):
        """ Get the efficiency using a logistic model for the 
            given event

            Args:
                sample: (PyWIMPs Sample)

            Returns: the efficiency
        """

        return self.eff/(1 +np.exp(-(sample.Er - self.x0)/self.xscale ) )
    
