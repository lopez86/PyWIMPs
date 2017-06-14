""" DetectorModel.py

    Holds a container class to save information about a detector's
    response to an event.

"""
__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

from .Efficiency import Efficiency
from .Response import Response
from ..mc.sample import Sample
import numpy as np

class DetectorModel:
    """ Container to hold detector response model 

        Attributes:
            efficiency: Detector efficiency model
            response: Detector response model
    """
    def __init__(self):
        """ Initializes with default attributes 
    
            The efficiency model is set to a constant.
            The response is set to do nothing.
            The random number generator is the default numpy one.
    
        """
        self.efficiency = Efficiency()
        self.response = Response()
        self._rand = np.random

    def initialize(self):
        """ Do any initial calculations of parameters."""
        self.efficiency.initialize() 
        self.response.initialize()


    def set_params(self,pars):
        """ Set parameters and pass to contained objects.
    
            This class has no parameters to be set and passes the
            parameters to the efficiency and response models. The
            random number generator must be set separately.
    
            Args:
               pars: (dict{string}) the new parameters
    
        """
        self.efficiency.set_params(pars)
        self.response.set_params(pars)

    @property
    def random(self):
        """ Numpy random number generator for this object"""
        return self._rand
    
    @random.setter
    def set_random(self,r):
        """ Sets the random number generator from a Numpy
            RandomState or equivalent.

            Args:
                r: (Numpy RandomState)
            
        """
        self._rand = r
        self.efficiency.set_random(r)
        self.response.set_random(r)
        


    def weighted_throw(self,sample):
        """ Applies the detector response and efficiency to
            a truth sample.

            Applies the efficiency and response. A detector
            weight is generated using the efficiency, while
            the response generates random reconstructed
            event properties.

            Args:
                sample: (PyWIMPs Sample)

            Returns:
                sample (PyWIMPs Sample)
        """
        sample.det_weight = self.efficiency.efficiency(sample) * sample.det_weight
        return self.response.weighted_throw(sample)

    def unweighted_throw(self,sample):
        """ Applies the detector response and efficiency to
            a truth sample.

            Applies the efficiency and response. A binary
            efficiency weight of 0 or 1 is generated and
            the detector response is thrown according to 
            the response function.

            Args:
                sample: (PyWIMPs Sample)

            Returns:
                sample (PyWIMPs Sample)
        """
        eff = self.efficiency.efficiency(sample)
        rnd = self._rand.rand()
        if eff < rnd:
            sample.det_weight = 0
        return self.response.unweighted_throw(sample)
