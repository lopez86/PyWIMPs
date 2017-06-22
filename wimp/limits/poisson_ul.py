""" poisson_ul.py

Module for simple Poisson-distribution based limits.
These assume some confidence level, and possibly an
expected number of background events. With no backgrounds,
the user will get a model-independent result.

Implemented here are:

* Frequentist Poisson upper limit assuming no backgrounds.
* CL_s upper limit for frequentist limits with a known background.
* Feldman-Cousins limit/interval to ensure coverage while 
  transitioning from an upper limit to an interval.
* Bayesian upper limit with a background model with a uniform
  prior. This is equivalent to the CL_s limit.
* Bayesian upper limit with a Jeffreys prior and background model.

The first 4 of these should be working correctly and some
validation is probably needed for the last.

TODO:
  Many of the methods here are also included in the other
  modules of this sub-package. I've found that it is better to
  have this type of limit setting as separate functions rather 
  than classes that need to know about the specific model.

"""

def freq_bkg_free_ul(N_exp=0,CL=0.9,tol=1e-4):
    """ Function to calculate Poisson upper limits.

        Args:
            N_exp: Number of measured events.
            CL: Confidence level
    """
    # We measured nothing. Easy case
    if N_exp == 0:
        # Works for frequentist and also 
        # Bayesian with flat prior
        # The limit in number of events
        n_limit = -np.log(1 - CL) 
        return n_limit

    # We actually measured some events
    err = 2
    fn_last = np.zeros(3) + 2
    n_limit = 0.5*np.sqrt(N_exp)+N_exp # First guess
    fn = 0
    while err > tol:
        fn = ( CL - 1)
        dfdn = 0
        for i in range(N_exp+1): # Should the +1 be here?
            pmf = scipy.stats.poisson.pmf(i,n_limit)
            fn = fn + pmf 
            dfdn = dfdn + (-1 + i / n_limit) * pmf
        ## Newton's method
        

        #dfdn = np.sign(dfdn) * np.max([np.abs(dfdn),0.01])
        n_limit = n_limit - fn / dfdn
        if n_limit <= 0:
            print('Error: Mean has fallen below 0')
            print(n_limit,fn,dfdn,fn/dfdn)
            return 0
        err = np.max( np.abs(fn )) / (1-CL)
        if np.abs( 1 - fn_last[2]/fn ) < (tol)**2 :
            n_limit = n_limit + tol
        fn_last[0] = fn_last[1]
        fn_last[1] = fn_last[2]
        fn_last[2] = fn
        
    return 0,n_limit

def cls_ul(N,b,CL=0.9, tol = 1e-4):
    """ Calculates CLs limits based on a counting experiment.
        
    Binned CLs limits are also used in collider experiments
    but for dark matter, that's probably not necessary. Might
    be a nice addition at some point.

    Args:
        N: The number of measured events
        b: The expected number of backgrounds
        CL: The confidence level
        tol: The tolerance used for stopping

    Returns:
        (s_low,s_high): The lower and upper limits of the interval
                        s_low = 0 here
    """
    pvb = scipy.stats.poisson.cdf(N,b)
    
    s = max(5,N-b) # 5 is just arbitrary here
    #nmax = 20
    it = 0
    while(1):#it<nmax):
       it+=1
       pois = np.array([scipy.stats.poisson.pmf(i,s+b) 
                        for i in range(0,N+1)])
       sum1 = np.sum(pois)
       CLs = sum1 / pvb
 #      print(s,CLs,sum1,pvb)
       diff = 1-CL - CLs
       dfds = -np.sum( (np.arange(0,N+1)/(1.0*s+b) -1)*pois )/pvb

       if diff/dfds < s:
           s = s - diff/dfds
       else:
           s = s * 0.5

       if np.abs(diff/(1.0-CL)) < tol:
           break

    return 0, s

def bayes_unif_ul(N_exp=0,b=0,CL=0.9,tol=1e-4):
    """ Returns a Bayesian upper limit with a uniform 
        prior.

    Identical to the CLs limit.

    Args:
        N: The number of measured events
        b: The expected number of backgrounds
        CL: The confidence level
        tol: The tolerance used for stopping

    Returns:
        (s_low,s_high): The lower and upper limits of the interval
                        s_low = 0 here
    """
    return cls_limit(N_exp,b,CL, tol)

def bayes_jeffreys_ul(N_exp=0,b=0,CL=0.9,tol=1e-4):
    """ Returns a Bayesian upper limit with a Jeffreys 
        prior.
    
    Args:
        N: The number of measured events
        b: The expected number of backgrounds
        CL: The confidence level
        tol: The tolerance used for stopping

    Returns:
        (s_low,s_high): The lower and upper limits of the interval
                        s_low = 0 here
    """

    err = 2
    fn_last = np.zeros(3) + 2
    s = 0.5*np.sqrt(N_exp)+N_exp # First guess

    g1_2 = scipy.special.gamma(N_exp+0.5)

    denom = scipy.special.gammaincc(N_exp+0.5,b)
    num_bkg = scipy.special.gammainc(N_exp+0.5,b)
    while(1): 
          
        integrand = exp(-(s+b)) * (s+b)**(N_exp+0.5)
        num_sb = scipy.special.gammainc(N_exp+0.5,s+b)
        fs = CL - (num_sb - num_bkg) / denom
        dfds = integrand / (g1_2 * denom)
        # Newton's method:
        s = s - fn / dfdn
        if n_limit <= 0:
            print('Error: Mean has fallen below 0')
            print(n_limit,fn,dfdn,fn/dfdn)
            return 0
         err = np.max( np.abs(fn )) / (1-CL)
        if np.abs( 1 - fn_last[2]/fn ) < (tol)**2 :
            s = s + tol
        fn_last[0] = fn_last[1]
        fn_last[1] = fn_last[2]
 
    return 0,s


def feldman_cousins(N_exp=0,b=0,CL=0.9):
    """ Function to calculate Feldman-Cousins confidence
        intervals.

        Args:
            N_exp: (int)  Number of measured events. Default: 0
            b: (float) Expected number of backgrounds. Default: 0
            CL: (float, 0-1) Confidence level. 
    """
    
    def fc_limits(s,bkg=0,CL0=0.9):
        """ Function to calculate bounds on the number of measured
            events for a model given the Feldman-Cousins ordering
            principle.
    
            Args:
                s: The expected number of signal events
                bkg: The expected number of background events
                CL0: The confidence interval
    
            Returns:
                The CL bounds on the number of expected events
                in the experiment given a signal and background
                model.
        """
        # Get approx. # of sigma:
        sigma = scipy.special.erfinv(CL0) * np.sqrt(2)
        N_max = int(max(s+b + 4*sigma * np.sqrt(s+bkg),20) )
        # Should be enough for a reasonable data set
        while True:
            N_max = N_max * 2
            ratio = np.zeros(N_max)
            pval = np.zeros(N_max)
            for n in range(N_max):
                pval[n] = scipy.stats.poisson.pmf(n,s+bkg)
                pbest = 1e4
                if b > n:
                    pbest = scipy.stats.poisson.pmf(n,bkg)
              
                    # minus sign is because sorting is in ascending order
                else:
                    pbest = scipy.stats.poisson.pmf(n,n)
                ratio[n] = -pval[n] / pbest  
            sorted_indices = np.argsort(ratio)
            total_prob = 0
            lim_min = N_max
            lim_max = -1
            n = 0
            while total_prob < CL0 and n < len(ratio):
                total_prob += pval[sorted_indices[n]]
                if sorted_indices[n] < lim_min:
                    lim_min = sorted_indices[n]
                if sorted_indices[n] > lim_max:
                    lim_max = sorted_indices[n]
                n = n+1
        
            if lim_max < N_max-1:
                break
      
        return lim_min,lim_max,total_prob

    low_limit = 0
    sigma = scipy.special.erfinv(CL) * np.sqrt(2)
    up_limit = int(max(N_exp+10*sigma*np.sqrt(N_exp),20))

    lim_minl,lim_maxl,total_probl = fc_limits(low_limit,b,CL)

    lim_minu,lim_maxu,total_probu = fc_limits(up_limit,b,CL)
    while lim_minu <= N_exp:
        up_limit = 2 * up_limit
        lim_minu,lim_maxu,total_probu = \
                 fc_limits(up_limit,b,CL)

    upper_lim = False
    
    if lim_minl <= N_exp < lim_maxl:
        upper_lim = True # We want an interval in this case
        print("Calculating upper limit")
    # Binary search for the limit
    if upper_lim == False:
        low_l = low_limit
        up_l = up_limit
        
        while up_l - low_l > FeldmanCousinsLimit._tol:
            new_l = 0.5 * (up_l - low_l) + low_l
            minl, maxl, total_prob = fc_limits(new_l,b,CL)
            if maxl >= N_exp:
                up_l = new_l
            else:
                low_l = new_l
        low_limit = low_l

    # Upper limit now

    low_l = 0
    up_l = up_limit
    
    while up_l - low_l > FeldmanCousinsLimit._tol:
        new_l = 0.5 * (up_l - low_l) + low_l
        minl, maxl, total_prob = fc_limits(new_l,Nb,CL)
        if minl <= N_exp:
            low_l = new_l
        else:
            up_l = new_l
    up_limit = low_l


 
