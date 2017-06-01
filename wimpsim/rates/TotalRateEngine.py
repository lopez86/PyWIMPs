from wimpsim.units import speed_of_light
from wimpsim.units import MeV
from wimpsim.units import keV
from wimpsim.units import km_s

class TotalRateEngine:

    def __init__(self):
      self.calc_err = False
      self.max_iter = 100000

    def calculate_impweight(self,int_model,astro_model):
        vE = astro_model.vE()
        v0 = astro_model.v0()
        vesc = astro_model.vesc()
        dens = astro_model.wimp_density()
        Mx = int_model.Mx()
        Mt = int_model.Mt()
        Mtot = int_model.Mtot()
        cs = int_model.total_xs()
        mu = Mx*Mt/(Mx+Mt)

        while niter < self.max_iter:

            vec = np.random.normal(-vE,v0/np.sqrt(2),3)
            vec2 = vec + vE

            #Probability to throw this from the full Maxwellian distribution
            vec_prob = 1./(np.pi*v0*v0)**1.5 * \
                       np.exp( - (vec2.dot(vec2)) / (v0*v0) ) 
            vec_norm = vec.dot(vec)
            Ex = 0.5 * Mx * vec_norm / (speed_of_light*speed_of_light)  
            Emax = int_model.cross_section().Emax(Ex)

            Er = np.random.rand() * Emax
            vec_prob = vec_prob / Emax #Add the probability for this Er throw

            # We are calculating the value of f(v)|F|^2/v

            fval = astro_model.velocity(vec) * \
                   int_model.form_factor(2*Mt*Er) \
                   /vec_norm / vec_prob
            

            fave = fave + fval

            niter = niter + 1
        
        fave = fave / niter
        A = cs * dens * Mtot / (2*mu*mu*Mx)
        return A * fave

    def calculate_sample(int_model,astro_model)
        vE = astro_model.vE()
        v0 = astro_model.v0()
        vesc = astro_model.vesc()
        dens = astro_model.wimp_density()
        Mx = int_model.Mx()
        Mt = int_model.Mt()
        Mtot = int_model.Mtot()
        cs = int_model.total_xs()
        mu = Mx*Mt/(Mx+Mt)



          
          
      
