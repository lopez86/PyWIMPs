""" Experiment.py

    Class to represent an entire experiment.

    Includes a detector mode, astrophysics model,
    cross section model, and sampling models.


"""
__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"
from .DetectorModel import DetectorModel
from ..astro.AstroModel import AstroModel
from ..xsec.InteractionModel import InteractionModel
from ..mc.MaxwellWeightedSampler import MaxwellWeightedSampler
from ..mc.UniformWeightedSampler import UniformWeightedSampler
from ..mc.AcceptRejectSampler import AcceptRejectSampler
from ..mc.sample import Sample
from .. import units
import numpy as np

class Experiment:
    """ Class to hold all the information for a single experiment.

        Includes the astrophysics, interaction, and detector models,
        weighted and unweighted sampling.

        Attributes:
            rate_sampler: A sampling class
            event_sampler: An unweighted sampling class
            exposure: Total run time
            Emin: Minimum energy to consider in rates and samples
            Emax: Maximum energy to consider in rates and samples
            Nsamples: Number of samples to use in rate calculations
            nrec_meas: Number of measured events in a dataset
            nrec_true: Number of true events in a dataset
            nrec_total: Number of true events in a dataset,   
                        regardless of any analysis cuts

            nrec_meas_err: Approximate error on nrec_meas
            nrec_true_err: Approximate error on nrec_true
            nrec_total_err: Approximate error on nrec_total
    """
    def __init__(self):
        """ Initialize to some default values. """
        self._rand = np.random
        self._astro_model = AstroModel()
        self._interaction = InteractionModel()
        self.rate_sampler = \
                 MaxwellWeightedSampler(self.astro_model, \
                                        self.interaction)
        self.event_sampler = \
                 AcceptRejectSampler(self.astro_model, \
                                     self.interaction)

        self._detector_model = DetectorModel()

        self.exposure = 300. * units.day
        self.Emin = 0
        self.Emax = 100 * units.keV
        self.Nsamples = 100000
        self.nrec_meas = 0
        self.nrec_true = 0
        self.nrec_total = 0

        self.nrec_meas_err = 0
        self.nrec_true_err = 0
        self.nrec_total_err = 0

    def set_params(self,pars):
        """ Set parameters from a dictionary.

            Args:
                pars: {string}

            Parameters:
                Exposure: total integrated run time
                ExpEmin: Minimum energy in analysis
                ExpEmax: Maximum energy in analysis
                ExpNSamples: Number of samples

        """
        if 'Exposure' in pars:
            self.exposure = pars['Exposure']
        if 'ExpEmin' in pars:
            self.Emin = pars['ExpEmin']
        if 'ExpEmax' in pars:
            self.Emax = pars['ExpEmax']
        if 'ExpNSamples' in pars:
            self.Nsamples = pars['ExpNSamples']
        self.detector_model.set_params(pars)
        self.astro_model.set_params(pars)
        self.interaction.set_params(pars)
        self.rate_sampler.set_params(pars)
        self.event_sampler.set_params(pars)

    def initialize(self):
        """ Do any initial calculations of parameters."""
        self.detector_model.initialize()
        self.astro_model.initialize()
        self.interaction.initialize()
        self.rate_sampler.initialize()
        self.event_sampler.initialize()

    @property
    def random(self):
        """ The random number generator. """
        return self._rand
 
    @random.setter
    def random(self,r):
        """ Set the random number generator.

            Args:
                r (Numpy RandomState or equivalent)
        """
        self._rand = r
        self.detector_model.random = r
        self.astro_model.random = r
        self.interaction.random = r
        self.rate_sampler.random = r
        self.event_sampler.random = r

    @property
    def astro_model(self):
        """ The astrophysics model. """
        return self._astro_model

    @astro_model.setter
    def set_astro(self,am):
        """ Set the astrophysics model. 

            Args:
                am (AstroModel)
        """
        self._astro_model = am 
        self.rate_sampler.astro_model = am
        self.event_sampler.astro_model = am

    @property
    def interaction(self):
        """ The interaction model."""
        return self._interaction

    @interaction.setter
    def set_interaction(self,im):
        """ Set the interaction model.

            Args:
               im (InteractionModel)
        """
        self._interaction = im
        self.rate_sampler.interaction = im
        self.event_sampler.interaction = im

    @property
    def detector_model(self):
        """ The detector model."""
        return self._detector_model

    @detector_model.setter
    def set_detector_model(self,det):
        """ Set the detector mode.

            Args:
                det (DetectorModel)
        """
        self._detector_model = det

    
    def initialize(self):
        """ Initialize the various pieces of the experiment."""
        self.rate_sampler.initialize()
        self.event_sampler.initialize()


    def event_rates(self,N = -1):
        """ Get the total, experimental, and true (in
            experimental bounds) rates from some number
            of throws of the rate sampler.

            Args:
                N (int): the number of throws. If 
                    not positive, use self.Nsamples

            Returns:
                A dictionary with the different weights.
                The keys are:
                Total, TotalErr, Truth, TruthErr, 
                Meas, and MeasErr

        """
        self.nrec_total = 0
        self.nrec_true = 0
        self.nrec_meas = 0

        self.nrec_total_err = 0
        self.nrec_true_err = 0
        self.nrec_meas_err = 0
        if N <= 0: 
            N = self.Nsamples

        for i in range(N):
            s = self.rate_sampler.sample()
            self.nrec_total += s.gen_weight
            self.nrec_total_err += s.gen_weight**2
            if self.Emin <= s.Er < self.Emax:
                self.nrec_true += s.gen_weight
                self.nrec_true_err += s.gen_weight**2 
            s = self.detector_model.weighted_throw(s)
            if self.Emin <= s.Er_reco < self.Emax:
                self.nrec_meas += s.weight 
                self.nrec_meas_err += s.weight**2

        #print("Integral: ",self.nrec_total)      
        self.nrec_total *= self.exposure / N
        self.nrec_true *= self.exposure / N
        self.nrec_meas *= self.exposure / N

        self.nrec_total_err = self.exposure * np.sqrt(self.nrec_total_err) /N
        self.nrec_true_err = self.exposure * np.sqrt(self.nrec_true_err) /N
        self.nrec_meas_err = self.exposure * np.sqrt(self.nrec_meas_err) /N
 

        return {'Total':self.nrec_total,
                'TotalErr':self.nrec_total_err,
                'Truth':self.nrec_true,
                'TruthErr':self.nrec_true_err,
                'Meas':self.nrec_meas,
                'MeasErr':self.nrec_meas_err}

    def throw_event(self):
        """ Throw a single unweighted event. 

            Returns:
               The sample for the event.
        """
        
        sample = None
        while True:
            sample = self.event_sampler.sample()
            sample = self.detector_model.unweighted_throw(s)
            if (s.weight > 0) and (self.Emin<=s.Er<self.Emax):
                break

        return sample

    def throw_dataset(self,aveN):
        """ Throw a whole dataset given an average
            number of events in the dataset.

            Returns:
                A list of events.
        """

        N = self.rand.poisson(aveN)
        return [data.throw_event() for i in range(N)]
