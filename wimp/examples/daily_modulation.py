""" daily_modulation.py

This example looks at the daily modulation of the signal
direction as Earth rotates around its axis. This gives 
a sidereal modulation (i.e. a full rotation rather than a
single 24-hour day), which can be distinguished from
backgrounds from the lab environment or from the sun.

Example:

>>> run_example()

This will make several plots.

""" 

__author__ = 'Jeremy P. Lopez'
__date__ = 'October 2017'
__copyright__ = '(c) October 2017, Jeremy P. Lopez'

from .. import units
from ..astro import Coordinates
from ..astro import locations
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap

def CygnusDirection(tmin,tmax,npts=48):
    """ Calculates the direction of the dark matter wind (roughly
        the same as the constellation Cygnus) for some set of points.
        The calculation is done for the SURF lab at the Homestake Mine
        in Lead, SD, USA.

    Args:
        tmin (Datetime): The start time
        tmax (Datetime): The end time
        npts: The number of points

    Returns: (alts, azs, times)
        alts: Altitudes in lab frame
        azs: Azimuths in lab frame
        times: Times as Unix timestamps

    """
    
    utc_start = dt.datetime(1970,1,1,0,0,0,tzinfo=dt.timezone.utc)
    start_time = (tmin-utc_start).total_seconds()
    end_time = (tmax - utc_start).total_seconds()
 
    pts = np.mgrid[start_time:end_time:complex(0,npts)]
    alts = np.zeros(npts)
    azs = np.zeros(npts)

    coord = Coordinates()
    loc = locations.SURF

    for i in range(npts):
        pt = pts[i]
        theTime = dt.datetime.fromtimestamp(pt,tz = dt.timezone.utc)
        vE = coord.earth_motion_lab(theTime,loc)
        uE = vE / np.sqrt(np.dot(vE,vE))
        azs[i] = 0.5 * np.pi - np.arctan2(uE[1],uE[0])
        alts[i] = 0.5 * np.pi - np.arccos(uE[2])

    return alts,azs,pts

def run_example():
    """ Runs the example and creates several plots of how the dark
        matter direction changes over time.
    """
    fig = plt.figure(1)
    fig.add_subplot(111)
    tmin = dt.datetime(2017,6,1,0,0,0,tzinfo=dt.timezone.utc)
    tmax = dt.datetime(2017,6,2,0,0,0,tzinfo=dt.timezone.utc)
    alts,azs,pts = CygnusDirection(tmin,tmax)
    pts =24*(pts - np.min(pts)) / (np.max(pts)-np.min(pts))
    bmap = Basemap(projection='hammer',lon_0=0)
    bmap.drawmeridians(np.arange(0,360,60),color='0.5')
    bmap.drawparallels(np.arange(-90,90,30),color='0.5')
    alts = 180./np.pi * alts
    azs = 180./np.pi * azs
    x,y = bmap(azs,alts)
    #print('Min altitude: ' + str(np.min(alts)))
    bmap.scatter(x,y,c=pts,marker='.',cmap='jet')
    cb = bmap.colorbar(location='right',label='Time of Day [hr]')
    plt.title('Source Direction at SURF on June 1st, 2017')

    cb.set_ticks([0,6,12,18,24.0])

    fig = plt.figure(2)
    fig.add_subplot(111)
    tmin = dt.datetime(2017,9,1,0,0,0,tzinfo=dt.timezone.utc)
    tmax = dt.datetime(2017,9,2,0,0,0,tzinfo=dt.timezone.utc)
    alts,azs,pts = CygnusDirection(tmin,tmax)
    pts =24*(pts - np.min(pts)) / (np.max(pts)-np.min(pts))
    bmap = Basemap(projection='hammer',lon_0=0)
    bmap.drawmeridians(np.arange(0,360,60),color='0.5')
    bmap.drawparallels(np.arange(-90,90,30),color='0.5')
    alts = 180./np.pi * alts
    azs = 180./np.pi * azs
    x,y = bmap(azs,alts)
    
    #print('Min altitude: ' + str(np.min(alts)))
    bmap.scatter(x,y,c=pts,marker='.',cmap='jet')
    cb = bmap.colorbar(location='right',label='Time of Day [hr]')
    plt.title('Source Direction at SURF on Jan. 1st, 2017')

    cb.set_ticks([0,6,12,18,24.0])
    plt.show()

