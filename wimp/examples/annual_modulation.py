""" This example calculates the annual modulation of SI recoil events for 50 GeV WIMPs on xenon """

from .. import units
from .. import xsec
from .. import astro
from .. import mc
from .. import det

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def AnnualModulation(tmin,tmax,npts):

    pars = {'AtomicNumber':131,\
            'Mx':50 * units.GeV,\
            'Mt':131 * units.amu,\
            'XS':2e-39 * units.cm*units.cm,\
            'EffMax':1.0,\
            'Exposure':1.0 * units.day,\
            'Mtot':100*units.kg,\
            'v0':230*units.km/units.sec,\
            'vesc':544*units.km/units.sec,\
            'rhox':0.3*units.GeV/(units.cm**3)
           }

    ex = det.Experiment()
    ex.set_params(pars)
    utc_start = dt.datetime(1970,1,1,0,0,0,tzinfo = dt.timezone.utc)
    start_time = (tmin - utc_start).total_seconds()
    end_time = (tmax - utc_start).total_seconds()

    pts = np.mgrid[start_time:end_time:complex(0,npts)]
    rates = np.zeros(npts)
    rate_errs = np.zeros(npts)
    v = np.zeros(npts)
    
    coord = astro.Coordinates()
    loc = astro.locations.SURF

    for i in range(npts):
        pt = pts[i]
        theTime = dt.datetime.fromtimestamp(pt,tz=dt.timezone.utc)
        vE = coord.earth_motion_lab(theTime,loc)
        v[i] = np.sqrt(vE.dot(vE)) * units.sec / units.km

        ex.set_params({'vE':vE})
        # We've changed the model so need to reinitialize
        ex.initialize()
        # Calculate the rates (here, events/day)
        rate = ex.event_rates()
        rates[i] = rate['Total'] 
        rate_errs[i] = rate['TotalErr'] 
        print(i,rates[i])

    return pts,rates,rate_errs,v

def run_example():
    tmin = dt.datetime(2016,1,1,0,0,0,tzinfo = dt.timezone.utc)
    tmax = dt.datetime(2018,1,1,0,0,0,tzinfo = dt.timezone.utc)
   
    times,rates,rate_errs,v = AnnualModulation(tmin,tmax,1)
    return
    dates = [dt.datetime.fromtimestamp(ts,tz=dt.timezone.utc) for ts in times]
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.plot(dates,v)
    ax.set_ylabel('WIMP Velocity [km/s]') 
    ax.set_xlim(dates[0],dates[-1])

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()


    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    ax.plot(dates,rates)
    ax.set_ylabel('Interactions per Day') 
    ax.set_xlim(dates[0],dates[-1])

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()

    fig = plt.figure(3)
    ax = fig.add_subplot(111)

    ax.plot(dates,rate_errs/rates)
    ax.set_ylabel('Fractional Uncertainty') 
    ax.set_xlim(dates[0],dates[-1])

    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')

    ax.grid(True)
    fig.autofmt_xdate()



    plt.show()
