""" systematics.py

This example shows how you might start to study systematics
using one of the simple counting experiment limit methods.

The data here is meant to approximate the parameters of the
LUX 2014-2016 run, which as of 2017 has produced the 
best WIMP-nucleon spin-independent limit of any direct
detection experiment.

This takes many throws of the parameters of a detector
model and produces a histogram and some statistics about
the resulting distribution of 90% CL upper limits for
a 50 GeV WIMP.


Example:

>>> systematics()

  Produces a histogram of limits for different detector
  models and some associated statistics.

"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

from ..det.Experiment import Experiment
from ..xsec.HelmFormFactor import HelmFormFactor
from ..det.Efficiency import Efficiency
from ..limits.UpperLimitBkgFree import upper_limit
from ..det.DetectorModel import DetectorModel
from ..det.Response import GaussianResponse
from ..xsec.InteractionModel import InteractionModel

from ..xsec.SINormalization import SINormalization
from ..xsec.Nucleus import Nucleus
from .. import units
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from multiprocessing import Pool

class lux_efficiency(Efficiency):
    """Very approximate empirical formula to fit the rough
       shape of the efficiency curve in their most recent paper
       except for the drop-off in the efficiency above ~50 keVr"""
    def __init__(self):
        Efficiency.__init__(self)
        self.p = [4.3,2.6,2.7]
    def efficiency(self,s):
        Er = s.Er
        if Er < 1 * units.keV:
            return 0
        Er = Er/units.keV
        # Logistic with a fast turn-on above 0
        return (1 - np.exp(- (Er/self.p[2])**4)) / \
               (1 + np.exp(- (Er - self.p[0])/self.p[1]))
    
def syst_throw(seed):
    """ This produces a rate for a given throw.
       
        Args:
            seed: The random seed used for this throw
    """
    A = 131
    pars = {'AtomicNumber':A,
            'XS':1e-40 * units.cm**2,
            'Mt':131*units.amu,
            'vE':254*units.km/units.sec * np.array([0,0,1]),
            'v0':230 * units.km/units.sec,
            'vesc':544 * units.km/units.sec,
            'Mx':50 * units.GeV, # Roughly the optimal mass
            'Mtot':100 * units.kg,
            'rhox':0.3 * units.GeV / (units.cm**3),
            'ExpEmin':1.0 * units.keV,
            'ExpEmax':50.0 * units.keV,
            'Exposure':332 * units.day,
            'RespSigma':0
           }

    rnd = np.random.RandomState(seed)
    pars['RespMean'] = rnd.normal(1,0.1)
    exper = Experiment()
    print('Process: '+str(seed) + ' Mean: ' + str(pars['RespMean']))
    ff = HelmFormFactor()
    eff = lux_efficiency()

    exper.detector_model.efficiency = eff
    exper.detector_model.response = GaussianResponse()
    exper.interaction.form_factor = ff
    exper.random = rnd
    exper.set_params(pars)
    exper.initialize()
    rate = exper.event_rates()['Meas']
    print('Seed: ' + str(seed) + ' Rate: '+str(rate))
    return rate

    

def systematics(pool_size=4):
    """ Produces the histogram and systematics.

        Args:
            pool_size: The size of the multiprocessing
                       pool. Typically the # of CPUs or less
    """
    A = 131
    pars = {'AtomicNumber':A,
            'XS':1e-40 * units.cm**2,
            'Mt':131*units.amu,
            'vE':254*units.km/units.sec * np.array([0,0,1]),
            'v0':230 * units.km/units.sec,
            'vesc':544 * units.km/units.sec,
            'Mx':50 * units.GeV, # Roughly the optimal mass
            'Mtot':100 * units.kg,
            'rhox':0.3 * units.GeV / (units.cm**3),
            'ExpEmin':1.0 * units.keV,
            'ExpEmax':50.0 * units.keV,
            'Exposure':332 * units.day,
            'RespSigma':0
           }

    nucl = Nucleus({'NuclMassNumber':131,
                    'NuclAtomicNumber':54,
                    'NuclMass':131*units.amu})
    sinorm = SINormalization()

    n_syst_throws = 500
    pool = Pool(pool_size)
    seeds = np.random.randint(0,1000000,n_syst_throws)
    result = pool.map(syst_throw,seeds)
    norm = sinorm.normalize(nucl,50 * units.GeV)
    xs = pars['XS']
    N_exp = 0
    limit = upper_limit(N_exp,0.9)
    xs_lim = np.array([xs * limit * norm / r for r in result])
#    print(xs_lim/units.cm**2)


# Result from 500 throws with a random starting seed to generate
# throw seeds
#    xs_lim = np.array([  
#   1.08411582e-46,   1.07992927e-46,   1.07738277e-46,   1.07649946e-46,
#   1.09920684e-46,   1.08602148e-46,   1.08503798e-46,   1.08841389e-46,
#   1.08980565e-46,   1.08283468e-46,   1.07983225e-46,   1.08318003e-46,
#   1.08895397e-46,   1.08107878e-46,   1.08354493e-46,   1.08706474e-46,
#   1.08252096e-46,   1.08916321e-46,   1.08173224e-46,   1.09284206e-46,
#   1.09603393e-46,   1.08905492e-46,   1.08567224e-46,   1.08068872e-46,
#   1.07503529e-46,   1.08075702e-46,   1.08099780e-46,   1.08247658e-46,
#   1.07664987e-46,   1.09746242e-46,   1.07412438e-46,   1.09057510e-46,
#   1.08222335e-46,   1.08620460e-46,   1.08356729e-46,   1.08523047e-46,
#   1.07927612e-46,   1.08384062e-46,   1.07940615e-46,   1.08265823e-46,
#   1.08293689e-46,   1.08817433e-46,   1.09781132e-46,   1.08471216e-46,
#   1.08127451e-46,   1.08994527e-46,   1.08091132e-46,   1.08139914e-46,
#   1.08403000e-46,   1.07765708e-46,   1.08232630e-46,   1.08034485e-46,
#   1.08088927e-46,   1.09467408e-46,   1.08250161e-46,   1.08165491e-46,
#   1.08358796e-46,   1.07933830e-46,   1.08248212e-46,   1.08588163e-46,
#   1.08507166e-46,   1.08632512e-46,   1.08186851e-46,   1.09901666e-46,
#   1.07958268e-46,   1.08310196e-46,   1.08136829e-46,   1.08568096e-46,
#   1.08305096e-46,   1.08522125e-46,   1.09001170e-46,   1.08457841e-46,
#   1.08449519e-46,   1.08507207e-46,   1.09745893e-46,   1.07953733e-46,
#   1.07802444e-46,   1.07976347e-46,   1.08935111e-46,   1.08199685e-46,
#   1.09370209e-46,   1.09142874e-46,   1.09067037e-46,   1.07930515e-46,
#   1.08268015e-46,   1.07867131e-46,   1.07832941e-46,   1.08201104e-46,
#   1.08321089e-46,   1.08688291e-46,   1.08654838e-46,   1.08214169e-46,
#   1.08524296e-46,   1.08187980e-46,   1.08268277e-46,   1.08214018e-46,
#   1.08012695e-46,   1.09005512e-46,   1.07447870e-46,   1.08299252e-46,
#   1.08034620e-46,   1.08116701e-46,   1.09228959e-46,   1.07291345e-46,
#   1.09274885e-46,   1.08320752e-46,   1.08264451e-46,   1.09226565e-46,
#   1.07892050e-46,   1.08251628e-46,   1.08544275e-46,   1.08038128e-46,
#   1.08620688e-46,   1.09795584e-46,   1.09408497e-46,   1.08213009e-46,
#   1.07972203e-46,   1.07847973e-46,   1.08851807e-46,   1.07709342e-46,
#   1.09292245e-46,   1.07734234e-46,   1.08252916e-46,   1.08315740e-46,
#   1.08650912e-46,   1.08048846e-46,   1.08057589e-46,   1.08838680e-46,
#   1.08011563e-46,   1.08709227e-46,   1.08275955e-46,   1.08278204e-46,
#   1.07938423e-46,   1.07963789e-46,   1.07854021e-46,   1.08159853e-46,
#   1.07591241e-46,   1.08387158e-46,   1.08023186e-46,   1.08429296e-46,
#   1.08366611e-46,   1.08907786e-46,   1.08708938e-46,   1.07766079e-46,
#   1.07651375e-46,   1.08427233e-46,   1.08778668e-46,   1.08210894e-46,
#   1.08252922e-46,   1.08186845e-46,   1.08684613e-46,   1.08445609e-46,
#   1.09139209e-46,   1.08363758e-46,   1.07786132e-46,   1.08984172e-46,
#   1.08499164e-46,   1.08912100e-46,   1.08699402e-46,   1.08237494e-46,
#   1.08400583e-46,   1.08323385e-46,   1.08386425e-46,   1.07663801e-46,
#   1.08718973e-46,   1.07844051e-46,   1.08149881e-46,   1.08328634e-46,
#   1.08104660e-46,   1.08829192e-46,   1.08422692e-46,   1.08559889e-46,
#   1.07938969e-46,   1.08229765e-46,   1.08028274e-46,   1.08485515e-46,
#   1.07907165e-46,   1.07626882e-46,   1.08541694e-46,   1.09330163e-46,
#   1.07934265e-46,   1.07938591e-46,   1.09559922e-46,   1.08643474e-46,
#   1.08206533e-46,   1.08184315e-46,   1.08398447e-46,   1.08427563e-46,
#   1.08099231e-46,   1.08101256e-46,   1.08141490e-46,   1.07937535e-46,
#   1.08524232e-46,   1.08625861e-46,   1.08228495e-46,   1.09148453e-46,
#   1.09349583e-46,   1.08005169e-46,   1.07974944e-46,   1.08648958e-46,
#   1.09404455e-46,   1.09013547e-46,   1.08936913e-46,   1.08864955e-46,
#   1.09118003e-46,   1.08214164e-46,   1.08543888e-46,   1.08283325e-46,
#   1.07952124e-46,   1.08284297e-46,   1.09257098e-46,   1.08550430e-46,
#   1.09104232e-46,   1.08169166e-46,   1.08166349e-46,   1.08424956e-46,
#   1.08476548e-46,   1.08935482e-46,   1.08236936e-46,   1.09284944e-46,
#   1.07940021e-46,   1.09059857e-46,   1.07720118e-46,   1.08440216e-46,
#   1.07937196e-46,   1.08076782e-46,   1.08025304e-46,   1.08517855e-46,
#   1.09410232e-46,   1.07765258e-46,   1.09342062e-46,   1.07976179e-46,
#   1.08742045e-46,   1.08161025e-46,   1.10628701e-46,   1.08587799e-46,
#   1.08999551e-46,   1.08666051e-46,   1.08742270e-46,   1.09296925e-46,
#   1.07958780e-46,   1.08526792e-46,   1.08643456e-46,   1.08491376e-46,
#   1.08217267e-46,   1.08191759e-46,   1.08071745e-46,   1.08043302e-46,
#   1.09719181e-46,   1.07813067e-46,   1.08468906e-46,   1.08092376e-46,
#   1.08331155e-46,   1.07793355e-46,   1.08310948e-46,   1.07680515e-46,
#   1.08581722e-46,   1.08469573e-46,   1.08028140e-46,   1.07732468e-46,
#   1.08790101e-46,   1.08161708e-46,   1.09183635e-46,   1.08167654e-46,
#   1.08421719e-46,   1.08678253e-46,   1.08033552e-46,   1.08000095e-46,
#   1.08569459e-46,   1.08301239e-46,   1.08853321e-46,   1.08190131e-46,
#   1.08216220e-46,   1.09237870e-46,   1.08128163e-46,   1.07925294e-46,
#   1.07626457e-46,   1.08178125e-46,   1.08382513e-46,   1.07946970e-46,
#   1.08271126e-46,   1.08119600e-46,   1.08399139e-46,   1.08167375e-46,
#   1.07981332e-46,   1.08242045e-46,   1.08825650e-46,   1.08650229e-46,
#   1.08121714e-46,   1.07813631e-46,   1.09772357e-46,   1.08311663e-46,
#   1.08972472e-46,   1.08544200e-46,   1.08621772e-46,   1.07893067e-46,
#   1.08609432e-46,   1.08286427e-46,   1.08478852e-46,   1.07884624e-46,
#   1.08028624e-46,   1.07416148e-46,   1.07848689e-46,   1.08495502e-46,
#   1.08164686e-46,   1.08355974e-46,   1.07908547e-46,   1.08771085e-46,
#   1.07949045e-46,   1.07509715e-46,   1.08060668e-46,   1.09257298e-46,
#   1.08250329e-46,   1.07906901e-46,   1.07707724e-46,   1.08037608e-46,
#   1.08641192e-46,   1.08520761e-46,   1.08996936e-46,   1.09079550e-46,
#   1.08015814e-46,   1.08077349e-46,   1.08372001e-46,   1.08105140e-46,
#   1.08224325e-46,   1.07765906e-46,   1.10677489e-46,   1.08135640e-46,
#   1.07942491e-46,   1.08269514e-46,   1.08682260e-46,   1.08208952e-46,
#   1.09010750e-46,   1.08016074e-46,   1.07721229e-46,   1.07133688e-46,
#   1.08152704e-46,   1.08245561e-46,   1.08444415e-46,   1.08360794e-46,
#   1.08461377e-46,   1.07926042e-46,   1.09233628e-46,   1.09296797e-46,
#   1.08400779e-46,   1.08349191e-46,   1.10529869e-46,   1.08361224e-46,
#   1.08802277e-46,   1.07918234e-46,   1.08111445e-46,   1.07970806e-46,
#   1.08110609e-46,   1.08241155e-46,   1.08313183e-46,   1.10016826e-46,
#   1.08599705e-46,   1.08487641e-46,   1.08924184e-46,   1.08145753e-46,
#   1.08011438e-46,   1.08062795e-46,   1.07879551e-46,   1.08996833e-46,
#   1.08233783e-46,   1.07841899e-46,   1.07506288e-46,   1.08763209e-46,
#   1.08537924e-46,   1.08034272e-46,   1.07961943e-46,   1.08880924e-46,
#   1.08185279e-46,   1.08939380e-46,   1.09023653e-46,   1.07842931e-46,
#   1.08055809e-46,   1.08136681e-46,   1.07858587e-46,   1.09113346e-46,
#   1.08608605e-46,   1.09032366e-46,   1.08353646e-46,   1.08423451e-46,
#   1.08662977e-46,   1.08199818e-46,   1.09502637e-46,   1.08236102e-46,
#   1.08938082e-46,   1.08172213e-46,   1.08557358e-46,   1.08151831e-46,
#   1.08097405e-46,   1.07824655e-46,   1.08009337e-46,   1.08173552e-46,
#   1.08837891e-46,   1.08207717e-46,   1.08132909e-46,   1.09385813e-46,
#   1.07985213e-46,   1.08777289e-46,   1.07969669e-46,   1.08435252e-46,
#   1.08396981e-46,   1.09555132e-46,   1.08023419e-46,   1.08949918e-46,
#   1.08405608e-46,   1.09538705e-46,   1.08844167e-46,   1.08272039e-46,
#   1.08629906e-46,   1.09456856e-46,   1.07855346e-46,   1.07522671e-46,
#   1.08659099e-46,   1.08571106e-46,   1.08404395e-46,   1.08454940e-46,
#   1.08450223e-46,   1.09853927e-46,   1.09344385e-46,   1.08796724e-46,
#   1.07783487e-46,   1.08102580e-46,   1.07696323e-46,   1.07896918e-46,
#   1.09778141e-46,   1.08085474e-46,   1.08272100e-46,   1.08500989e-46,
#   1.08989125e-46,   1.09531868e-46,   1.09227315e-46,   1.09284149e-46,
#   1.08818580e-46,   1.09003455e-46,   1.07971592e-46,   1.08162935e-46,
#   1.08845092e-46,   1.07913170e-46,   1.09595389e-46,   1.08562957e-46,
#   1.07943963e-46,   1.09659001e-46,   1.08018942e-46,   1.08242053e-46,
#   1.07690073e-46,   1.08006616e-46,   1.08042255e-46,   1.08360226e-46,
#   1.07995075e-46,   1.08269265e-46,   1.08322072e-46,   1.08048460e-46,
#   1.07807188e-46,   1.08707553e-46,   1.08852217e-46,   1.08287946e-46,
#   1.08152804e-46,   1.08362289e-46,   1.08444058e-46,   1.07547684e-46,
#   1.07918307e-46,   1.08767512e-46,   1.07731364e-46,   1.08916047e-46,
#   1.08174176e-46,   1.08564131e-46,   1.09991126e-46,   1.09551828e-46,
#   1.08041890e-46,   1.08395302e-46,   1.08494806e-46,   1.08023297e-46,
#   1.08977322e-46,   1.08378754e-46,   1.08434658e-46,   1.08241954e-46,
#   1.08628545e-46,   1.09243507e-46,   1.08601659e-46,   1.09087259e-46,
#   1.08396219e-46,   1.08406510e-46,   1.09003819e-46,   1.08567154e-46,
#   1.08379662e-46,   1.08250097e-46,   1.07606176e-46,   1.07818660e-46,
#   1.07894498e-46,   1.08322210e-46,   1.08497619e-46,   1.08220789e-46,
#   1.08024118e-46,   1.08667485e-46,   1.08959861e-46,   1.08463059e-46])



    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.hist(xs_lim / units.cm**2,bins=20,range=(1.070e-46,1.11e-46))
    ave = np.mean(xs_lim) 
    med = np.median(xs_lim)
    std = np.std(xs_lim)
    perc = np.percentile(xs_lim,[2.5,18,82,97.5])
    print('Mean: %.4e'%(ave,) )
    print('Median: %.4e'%(med,))
    print('Standard Deviation: %.2e'%(std,))
    print('68%% Band: %.4e %.4e'%(perc[1],perc[2]))
    print('95%% Band: %.4e %.4e'%(perc[0],perc[3]))
    
    ax.set_xlim([1.070e-46,1.11e-46])
    ax.text(1.09e-46, 65, '50 GeV WIMP on Xenon\nLUX 2014-16 Exposure', fontsize=13)
    ax.set_xlabel(r'90% CL Upper Limit [cm$^2$]',fontsize=15)
    ax.set_ylabel(r'Number of Throws',fontsize=15)
    ax.set_title('Effect of 10% Energy Scale Uncertainty',fontsize=15)
    ax.grid(alpha=0.3)
    plt.show()
