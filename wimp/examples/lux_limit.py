from ..det.Experiment import Experiment
from ..xsec.HelmFormFactor import HelmFormFactor

class lux_efficiency(Efficiency):
    """Very approximate empirical formulat to fit the rough
       shape of the efficiency curve in their most recent paper
       except for the drop-off in the efficiency above ~50 keVr"""
    def __init__(self):
        self.p = [4.3,2.6,2.7]
    def get_eff(self,s):
        Er = s.Er
        if Er < 1 * units.keV:
            return 0
        # Logistic with a fast turn-on above 0
        return (1 - np.exp(- (Er/self.p[2])**4)) / \
               (1 + np.exp(- (Er - self.p[0])/self.p[1]))
    
# We'll just leave the response as a Gaussian for now. That
# shouldn't matter too much for the limits


def lux_limits():
    A = 131
    ff = HelmFormFactor()
    ff.set_params({'A':A})
    im = InteractionModel()
    im.form_factor = ff

    vel = VelocityDist()
    vel.vE = np.array([0,0,254*units.km/units.sec])
    # Direction doesn't matter here
    vel.vesc = 544 * units.km / units.sec
    vel.v0 = 230
    fid_mass = 100 * units.kg # approximate
    rho = 0.3 * units.GeV / units.cm**3
    xs0 = 1 * cm**2 # We'll get cross sections in cm^2




