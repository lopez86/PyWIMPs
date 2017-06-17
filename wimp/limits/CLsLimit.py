""" CLsLimit.py

Calculate CL_s limits for a simple counting experiment.
"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'

import scipy.stats
import numpy as np


def cls_limit(N,b,CL=0.9, tol = 1e-4):
    """ Calculates CLs limits based on a counting experiment.
        
    Binned CLs limits are also used in collider experiments
    but for dark matter, that's probably not necessary. Might
    be a nice addition at some point.
    """
    pvb = scipy.stats.poisson.cdf(N,b)
    
    s = max(5,N-b) # 5 is just arbitrary here
    nmax = 20
    it = 0
    while(it<nmax):
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
