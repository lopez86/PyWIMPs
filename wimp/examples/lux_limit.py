from ..det.Experiment import Experiment
from ..xsec.HelmFormFactor import HelmFormFactor
from ..det.Efficiency import Efficiency
from ..limits.UpperLimitBkgFree import upper_limit
from ..det.DetectorModel import DetectorModel
from ..xsec.InteractionModel import InteractionModel

from ..xsec.SINormalization import SINormalization
from ..xsec.Nucleus import Nucleus
from .. import units
import numpy as np

class lux_efficiency(Efficiency):
    """Very approximate empirical formulat to fit the rough
       shape of the efficiency curve in their most recent paper
       except for the drop-off in the efficiency above ~50 keVr"""
    def __init__(self):
        self.p = [4.3,2.6,2.7]
    def efficiency(self,s):
        Er = s.Er
        if Er < 1 * units.keV:
            return 0
        Er = Er/units.keV
        # Logistic with a fast turn-on above 0
        return (1 - np.exp(- (Er/self.p[2])**4)) / \
               (1 + np.exp(- (Er - self.p[0])/self.p[1]))
    
# We'll just leave the response as a Gaussian for now. That
# shouldn't matter too much for the limits


def lux_limits():
    A = 131
    pars = {'AtomicNumber':A,
            'XS':1e-40 * units.cm**2,
            'Mt':131*units.amu,
            'vE':254*units.km/units.sec * np.array([0,0,1]),
            'v0':230 * units.km/units.sec,
            'vesc':544 * units.km/units.sec,
            'Mtot':100 * units.kg,
            'rhox':0.3 * units.GeV / (units.cm**3),
            'ExpEmin':1.0 * units.keV,
            'ExpEmax':50.0 * units.keV,
            'Exposure':332 * units.day,
           }

    nucl = Nucleus({'NuclMassNumber':131,
                    'NuclAtomicNumber':54,
                    'NuclMass':131*units.amu})
    sinorm = SINormalization()

    exper = Experiment()

    ff = HelmFormFactor()
    eff = lux_efficiency()

    exper.detector_model.efficiency = eff
    exper.interaction.form_factor = ff

    exper.set_params(pars)


    mass_grid = np.exp(np.arange(1,4,0.1) * np.log(10)) * units.GeV

    # Background free:
    N_exp = 0
    limit = upper_limit(N_exp,0.9)
    xs = exper.interaction.total_xs
    for m in mass_grid:
        exper.set_params({'Mx':m})
        exper.initialize()        
        rate = exper.event_rates()['Meas']
        xs_lim = xs * limit/rate
        print('Mass: ' + str(m/units.GeV) + ' GeV, XS: ' + str(xs_lim/units.cm**2)+' cm^2')
        xs_lim_norm = xs_lim * sinorm.normalize(nucl,m)
        print('\tNormalized: ' + str(xs_lim_norm/units.cm**2) + ' cm^2')


