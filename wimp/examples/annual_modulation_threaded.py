""" annual_modulation_threaded.py

This example calculates the annual modulation of SI recoil events 
for 50 GeV WIMPs on xenon. No form factor or detector effects
are considered here.

This is identical to annual_modulation.py except it
uses the multiprocessing package to create multiple threads
for the calculations. This requires being careful with things
like random number generators, which by default are global 
objects and probably aren't thread-safe.

Should run much faster than the non-threaded version.

Example:

>>> run_example()

Matplotlib will display several plots of things like
the total rate over time and the velocity.
"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) June 2017, Jeremy P. Lopez'

from .. import units
from .. import xsec
from .. import astro
from .. import mc
from .. import det

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from multiprocessing import Pool

# Fixed seeds here to always get same results
# One seed per thread
seeds = [198201, 418776, 471693, 263639, 888106, 722914, 188371]
vpts = np.array([225,230,235,240,245,250,260]) * units.km/units.sec

def calc_rate(i):
    """ We get the rates by interpolating from the velocity. 
    
    This function calculates the rate from one of our pre-set
    velocity points. It also creates a new random number generator
    to ensure thread safety.

    Args:
        i: The index of the velocity point and random number seed

    Returns:
        The ineteraction rate
    """
    rnd = np.random.RandomState(seeds[i])
    ex = det.Experiment()
    # For thread safety
    ex.set_random(rnd)
    vE = np.array([0,0,vpts[i]])
    pars = {'AtomicNumber':131,\
            'Mx':50 * units.GeV,\
            'Mt':131 * units.amu,\
            'XS':2e-39 * units.cm*units.cm,\
            'EffMax':1.0,\
            'Exposure':1.0 * units.day,\
            'Mtot':100*units.kg,\
            'v0':230*units.km/units.sec,\
            'vesc':544*units.km/units.sec,\
            'rhox':0.3*units.GeV/(units.cm**3),\
            'ExpNSamples':100000,\
            'vE':vE
           }
    ex.set_params(pars)
    ex.initialize()
    # Calculate the rates (here, events/day)
    rate = ex.event_rates()

    print(i,rate['Total'],rate['TotalErr'])
    return i,rate['Total'],rate['TotalErr']

def AnnualModulation(tmin,tmax,npts,pool_size=4):
    """ Calculates the event rate for the given time period.

    Args:
        tmin (Datetime): The start time
        tmax (Datetime): The end time
        npts: The number of points
        pool_size: The number of processes to run concurrently
                   Default: 4

    Returns: (note: return values are pretty ugly here)
        Time points
        Interaction rates (per day)
        Estimated rate uncertainties
        WIMP velocities
        Velocity points used for interpolation
        Rates at the interpolation points
        Rate errors at the interpolation points
   
    """


    utc_start = dt.datetime(1970,1,1,0,0,0,tzinfo = dt.timezone.utc)
    start_time = (tmin - utc_start).total_seconds()
    end_time = (tmax - utc_start).total_seconds()

    pts = np.mgrid[start_time:end_time:complex(0,npts)]
    rates = np.zeros(npts)
    rate_err = np.zeros(npts)
    v = np.zeros(npts)
    
    coord = astro.Coordinates()
    loc = astro.locations.SURF


    rate_pts = np.zeros(len(vpts))
    rate_idx = np.arange(0.0,len(vpts))
    rate_pt_err = np.zeros(len(vpts))
    pool = Pool(pool_size)
    m = pool.map( calc_rate,range(len(vpts)))

    for x in m:
        rate_pts[x[0]] = x[1]
        rate_pt_err[x[0]] = x[2]

    for i in range(npts):
        pt = pts[i]
        theTime = dt.datetime.fromtimestamp(pt,tz=dt.timezone.utc)
        #vE = np.array([0,0,230*units.km/units.sec])
        vE = coord.earth_motion_lab(theTime,loc)
        v[i] = np.sqrt(vE.dot(vE)) * units.sec / units.km

        rates[i] = np.interp(v[i] * units.km/units.sec,vpts,rate_pts)
        interp_val = np.interp(v[i] * units.km/units.sec,vpts,rate_idx)
        idx0 = int(interp_val)
        frac = interp_val % 1
        rate_err[i] = np.sqrt( frac**2 * rate_pt_err[idx0+1]**2 
                               + (1-frac)**2 * rate_pt_err[idx0]**2)
        
    return pts,rates,rate_err,v,vpts,rate_pts,rate_pt_err

def run_example():
    """ Gets the annual modulation signal for 2016 and 2017.

    Runs AnnualModulation()
    Creates a number of matplotlib plots and displays them.
    """
    tmin = dt.datetime(2016,1,1,0,0,0,tzinfo = dt.timezone.utc)
    tmax = dt.datetime(2018,1,1,0,0,0,tzinfo = dt.timezone.utc)
   
    times,rates,rate_err,v,vpts,rate_pts,rate_pt_err = \
                AnnualModulation(tmin,tmax,200)
    #return
    dates = [dt.datetime.fromtimestamp(ts,tz=dt.timezone.utc) 
             for ts in times]
    vpts = vpts * units.sec / units.km
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.plot(dates,v)
    ax.set_ylabel('WIMP Velocity [km/s]') 
    ax.set_xlabel('Date')
    ax.set_xlim(dates[0],dates[-1])

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()


    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    ax.plot(dates,rates)
    ax.fill_between(dates,rates-rate_err,rates+rate_err,alpha=0.3)
    ax.set_ylabel('Interactions per Day') 
    ax.set_xlabel('Date')
    ax.set_xlim(dates[0],dates[-1])

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()

    fig = plt.figure(3)
    ax = fig.add_subplot(111)

    ax.plot(vpts,rate_pts)
    ax.fill_between(vpts,rate_pts-rate_pt_err,rate_pts+rate_pt_err,alpha=0.3)
    ax.set_ylabel('Interactions per day') 
    ax.set_xlabel('WIMP Velocity [km/s]')
    ax.set_xlim(vpts[0],vpts[-1])

#    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)

    fig = plt.figure(4)
    ax = fig.add_subplot(111)

    ax.plot(dates,rates)
    ax.fill_between(dates,rates-rate_err,rates+rate_err,alpha=0.3)
    ax.set_ylabel('Interactions per Day') 
    ax.set_xlabel('Date')
    ax.set_xlim(dates[0],dates[-1])
    ax.set_ylim(0,np.max(rates+rate_err))

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()

    fig = plt.figure(5)
    ax = fig.add_subplot(111)

    ax.plot(vpts,rate_pts)
    ax.fill_between(vpts,rate_pts-rate_pt_err,rate_pts+rate_pt_err,alpha=0.3)
    ax.set_ylabel('Interactions per day') 
    ax.set_xlabel('WIMP Velocity [km/s]')
    ax.set_xlim(vpts[0],vpts[-1])
    ax.set_ylim(0,np.max(rate_pts+rate_pt_err))

#    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()

    fig = plt.figure(6)
    ax = fig.add_subplot(111)

    ax.plot(dates,v)
    ax.set_ylabel('WIMP Velocity [km/s]') 
    ax.set_xlabel('Date')
    ax.set_xlim(dates[0],dates[-1])
    ax.set_ylim(0,np.max(v))

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()

    plt.show()
