from DetectorModel import DetectorModel
from ..astro.AstroModel import AstroModel
from ..xsec.InteractionModel import InteractionModel
from ..mc.MaxwellWeightedSampler import MaxwellWeightedSampler
from ..mc.AcceptRejectSample import AcceptRejectSampler
from ..mc.sample import Sample
import ..units

class Experiment:
    def __init__(self):
        self.astro_model = AstroModel()
        self.interaction = InteractionModel()
        self.rate_sampler = \
                 MaxwellWeightedSampler(self.astro_model, \
                                        self.interaction)
        self.event_sampler = \
                 AcceptRejectSampler(self.astro_model \
                                     self.interaction)

        self.detector_model = DetectorModel()

        self.exposure = 300. * day
        self.Emin = 0
        self.Emax = 100 * units.keV
        self.nrec_meas = 0
        self.nrec_true = 0
        self.nrec_total = 0

        self.nrec_meas_err = 0
        self.nrec_true_err = 0
        self.nrec_total_err = 0


    def set_astro(self,am):
        self.astro_model = am 
        self.rate_sampler.astro_model = am
        self.event_sampler.astro_model = am

    def set_interaction(self,im):
        self.interaction = im
        self.rate_sampler.interaction = im
        self.event_sampler.interaction = im

    def event_rates(self,N = 100000):
        for i in range(N):
            s = self.rate_sampler.sample()
            self.nrec_total += s.weight
            self.nrec_total_err += s.weight**2
            if self.Emin <= s.Er < self.Emax:
                self.nrec_true += s.weight
                self.nrec_true_err += s.weight**2
            s = self.detector.weighted_throw(s)
            if self.Emin <= s.Er < self.Emax:
                self.nrec_meas += s.weight
                self.nrec_meas_err += s.weight**2

        self.nrec_total *= self.exposure
        self.nrec_true *= self.exposure
        self.nrec_meas *= self.exposure

        self.nrec_total_err = self.exposure * np.sqrt(self.nrec_total_err)
        self.nrec_true_err = self.exposure * np.sqrt(self.nrec_true_err)
        self.nrec_meas_err = self.exposure * np.sqrt(self.nrec_meas_err)

        return {'Total':self.nrec_total,
                'TotalErr':self.nrec_total_err,
                'Truth':self.nrec_true,
                'TruthErr':self.nrec_true_err,
                'Meas':self.nrec_meas,
                'MeasErr':self.nrec_meas_err}

    def throw_dataset(self,aveN):
        data = []
        N = np.random.poisson(aveN)

        while len(data) < N:
            s = self.event_sampler.sample()
            s = self.detector.unweighted_throw(s)
            if (s.weight > 0) and (self.Emin<=s.Er<self.Emax):
                data.append(s)

        return data
