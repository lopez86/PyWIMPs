""" Efficiency.py

    Basic classes to represent the detector response.


"""
__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"


import numpy as np
from ..mc.sample import Sample

class Response:
    """ Base detector response class: no change. """

    def __init__(self):
        """ Initialize the object. """
        self._rand = np.random
 
    def initialize(self):
        """ Do any initial calculations of parameters."""
        pass


    @property
    def random(self):
        """ The random number generator. """
        return self._rand

    def set_random(self,r):
        """ Set the random number generator.

            Args:
                r (Numpy RandomState or equivalent)

        """
        self._rand = r

    def set_params(self,pars):
        """ Set any parameters from a dictionary.
         
            Args:
                pars: {string}
        """
        pass

    def unweighted_throw(self,sample):
        """ Perform random throw over the distribution and
            return a sample with detector effects added.
            Good for generating realistic datasets.
    
            Args:
                sample (PyWIMPs sample)
    
            Returns:
                sample with detector response added
    
        """
        return sample

    def weighted_throw(self,sample):
        """ Perform random throw over the parameter space and
            return a sample with detector effects added. Here,
            throws do not need to follow the distribution and
            thus need an extra weight.
            Good for generating rates and distributions.
    
            Args:
                sample (PyWIMPs sample)
    
            Returns:
                sample with detector response added
    
        """
        return sample


class GaussianResponse(Response):
    """ The reconstructed energy follows a Gaussian distribution 
        with a width proportional to the true energy.

        Args:
            sigma (float): the fractional width of the distribution

    """
    def __init__(self):
        """ Initialize with a 10% resolution."""
        self.random = np.random
        self.sigma = 0.1

    def set_params(self,pars):
        """ Set parameters from a dictionary .
    
            Args:
                pars {string}
    
            Parameters:
                RespSigma: The fractional energy uncertainty
        """
        if 'RespSigma' in pars.keys():
            self.sigma = pars['RespSigma']

    def unweighted_throw(self,s):
        """ Perform random throw over the distribution and
            return a sample with detector effects added.
            Good for generating realistic datasets.
    
            Args:
                sample (PyWIMPs sample)
    
            Returns:
                sample with detector response added
    
        """
        s.Er_reco = self._rand.normal(s.Er,self.sigma * s.Er)
        return s

    def weighted_throw(self,s):
        """ Perform random throw over the parameter space and
            return a sample with detector effects added. Here,
            throws do not need to follow the distribution and
            thus need an extra weight.
            Good for generating rates and distributions.
    
            Args:
                sample (PyWIMPs sample)
    
            Returns:
                sample with detector response added
    
        """
        s.Er_reco = self._rand.normal(s.Er,self.sigma * s.Er)
        return s
