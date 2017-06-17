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

import matplotlib as mpl
import matplotlib.pyplot as plt

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


    mass_grid = np.exp(np.linspace(1,4,31) * np.log(10)) * units.GeV

    xs_pts = np.array(np.zeros(mass_grid.size))
    # Background free:
    N_exp = 0
    limit = upper_limit(N_exp,0.9)
    xs = exper.interaction.total_xs
    for i in range(mass_grid.size):
        m = mass_grid[i]
        exper.set_params({'Mx':m})
        exper.initialize()        
        rate = exper.event_rates()['Meas']
        xs_lim = xs * limit/rate
        print('Mass: ' + str(m/units.GeV) + ' GeV, XS: ' + str(xs_lim/units.cm**2)+' cm^2')
        xs_lim_norm = xs_lim * sinorm.normalize(nucl,m)
        print('\tNormalized: ' + str(xs_lim_norm/units.cm**2) + ' cm^2')
        xs_pts[i] = xs_lim_norm
#
#    print(mass_grid)
#    print(xs_pts)


#    xs_pts = np.array([  1.24296194e-45,   5.47519831e-46,   2.92874654e-46,   1.88419866e-46,
#   1.40448798e-46,   1.17375261e-46,   1.08433872e-46,   1.08367388e-46,
#   1.16197258e-46,   1.30532407e-46,   1.52182663e-46,   1.79776221e-46,
#   2.19920437e-46,   2.68062477e-46,   3.32156842e-46,   4.10870835e-46,
#   5.06458735e-46,   6.35652782e-46,   8.02296518e-46,   9.95962830e-46,
#   1.26033390e-45,   1.59625772e-45,   1.97989473e-45,   2.49643466e-45,
#   3.11955720e-45,   3.89523191e-45,   4.96394399e-45,   6.23846099e-45,
#   7.90421820e-45,   9.91687816e-45,   1.24262264635e-44])

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.plot(mass_grid/units.GeV,xs_pts/ units.cm**2)
    ax.fill_between(mass_grid/units.GeV,xs_pts/ units.cm**2,10**-43,alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([10,10000])
    ax.set_ylim([10**-47,10**-43])
    ax.text(20, 2e-44, '90% CL Excluded Region', fontsize=13)
    ax.text(500, 5e-47, 'Allowed Region', fontsize=13)
    ax.set_xlabel('WIMP Mass [GeV]',fontsize=15)
    ax.set_ylabel(r'WIMP-proton SI Cross Section [cm$^2$]',fontsize=15)
    ax.set_title('Estimated Sensitivity of LUX 2014-2016 Run',fontsize=15)
    ax.grid(alpha=0.3)
    plt.show()
