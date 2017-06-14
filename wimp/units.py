""" units.py

    This module holds the definitions of many units and physical
    constants. In the unit system used here, 

    Mass, energy: MeV = 1
    Time:         second = 1
    Distance:     cm = 1
    Angle:        rad = 1

"""

__author__    = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

import numpy as np

MeV = 1.
GeV = 1000.
keV = 0.001
kg = 5.60958e29
gram = 0.001*kg
tonne = 1000*kg
eV = 1e-6
TeV = 1e6
amu = 931.5 * MeV
cm = 1
fm = 1e-13
meter = 100
km = 1e5
sec = 1
minute = 60 * sec
hour = 3600 * sec
day = 86400 * sec
week = 7 * day
year = 365.25 * day
speed_of_light = 299792458 * meter / sec
hbarc = 197.327 * MeV * fm  # MeV-fm
deg = np.pi / 180.0
