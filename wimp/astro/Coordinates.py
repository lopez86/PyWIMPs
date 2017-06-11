import numpy as np
import datetime as dt
import astropy.time as astime
import astropy.coordinates as ascoords
from astropy import units as u
from .. import units

class Coordinates:
    """ 
    Class to calculate things like the direction and mean velocity of WIMPs
    in the dark matter halo.

    Equations for motion follow:
    Lewin, J.D. and Smith, P.F., Astropart. Phys. 6 (1996) 87-112.
    This is still referenced in more recent sources such as:
    Green, A.M. Phys. Rev. D68 (2003) 023004.
    and
    Freese, Lisanti, and Savage. Rev. Mod. Phys. 85 (2013) 1561.


    """ 
    def __init__(self):
        self.ur = np.array([0,230,0]) * units.km /units.sec
        self.us = np.array([9,12,7]) * units.km / units.sec
        self.beta = np.array([-5.503,59.575,29.812]) * units.deg
        self.lambda_ = np.array([266.141,-13.3485,179.3212]) * units.deg
        self.uEmean = 29.79 * units.km / units.sec
        self.ellipticity = 0.016722
        self.lambda0 = 13 * units.deg
        self.lambdaL0 = 280.460 * units.deg
        self.lambdaL1 = 0.9856474 * units.deg
        self.lambdag0 = 357.528 * units.deg
        self.lambdag1 = 0.9856003 * units.deg
        self.lambdaB = 1.915 * units.deg
        self.lambdaC = 0.020 * units.deg
        self.J2000 = dt.datetime(1999,12,31,12,0,0,tzinfo = dt.timezone.utc)

        self.earthSpeed = 6371. * units.km / (86400. * units.sec)

    def earth_motion_gal(self,timestamp):
        """ 
        Get the motion of Earth in galactic coordinates 
        A datetime object is needed in order to get the position
        around the sun


        Args:
            timestamp (datetime): The timestamp at which we want to do
                the calculation

        Returns:
            NumPy array of length 3. The velocity of Earth in galactic
            coordinates.

        """
        timedelta = timestamp - self.J2000#dt.datetime.fromtimestamp(timestamp) - self.J2000
        days = timedelta.total_seconds() / 86400.
        L = self.lambdaL0 + self.lambdaL1 * days
        g = self.lambdag0 + self.lambdag1 * days
        lambda_t = L + self.lambdaB * np.sin(g) + self.lambdaC * np.sin(2*g)
        uE_l = self.uEmean * (1 - self.ellipticity * np.sin(lambda_t - self.lambda0 ) )
        uE = uE_l * np.cos(self.beta) * np.sin(lambda_t - self.lambda_)
        return self.ur + self.us + uE

    def earth_motion_lab(self,timestamp,coord):
        """
        Get the velocity of the lab frame in local coordinates. 
        The coordinate system here is:
        (r,theta,phi) = (1,0,x) = up
        (r,theta,phi) = (1,90 deg, 0) = east
        (r,theta,phi) = (1,90 deg, 90 deg) = north

        This uses AstroPy to convert the galactic frame velocity
        into the altitude/azimuth system for a given time and 
        position. It then adds a small correction from the 
        rotation of Earth. There will be a very slight 
        sidereal modulation of the velocity.

        Args:
            timestamp (datetime): The time at which we want to
                calculate the velocity
            coord (list or array): The coordinates on Earth at 
                which we want to do the calculation. Of the 
                form [lat,lon].

        Returns:
            NumPy array of length 3. The velocity of the lab frame
            at the given time and place.

        """
        lat,lon = coord
        loc = ascoords.EarthLocation(lon=lon*u.deg,lat=lat*u.deg)
        time = astime.Time(timestamp) # A datetime
        earth_gal = self.earth_motion_gal(timestamp)
        earth_vel = np.sqrt(earth_gal.dot(earth_gal))
        gal_lon = np.arctan2(earth_gal[1],earth_gal[0])
        gal_lat = 0.5 * np.pi - np.arccos(earth_gal[2] / earth_vel)


        gal = ascoords.Galactic(l = gal_lon * u.rad,b = gal_lat * u.rad)
        gal_0 = ascoords.SkyCoord(gal)
        altaz = ascoords.AltAz(obstime=time,location=loc,pressure=0)
        altaz_coords = gal_0.transform_to(altaz)
        
        # Let's let x = east, y = north, z = up

        az = altaz_coords.data.lon.to(u.rad).value
        alt = altaz_coords.data.lat.to(u.rad).value
        v_earth = np.zeros(3)
        theta = 0.5 * np.pi - alt
        phi = 0.5 * np.pi - az
        v_earth[0] = np.sin(theta) * np.cos(phi) * earth_vel
        v_earth[1] = np.sin(theta) * np.sin(phi) * earth_vel
        v_earth[2] = np.cos(theta) * earth_vel


        # The lab frame is also traveling west due to Earth's rotation
        # Earth's rotation ends up being:
        v_rot = np.array([-self.earthSpeed,0,0])
        v_tot = v_earth + v_rot
        return v_tot
