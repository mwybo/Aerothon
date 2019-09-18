#!/usr/bin/env python

"""
Copyright (c) 2006, Russell A. Paielli <http://RussP.us>
All rights reserved.

This python script defines a basic system of physical units for air
traffic management (ATM).
"""

from __future__ import division # let 5/2 = 2.5 rather than 2

from os import environ as _environ

if "scalar_off" in _environ: from scalar_off import *
else: from scalar import *

rad = unit() # radians (dimensionless)
deg = unit("deg", rad * pi/180)

sec = unit("sec") # second ("s" is too short!)
Min = unit("Min", 60 * sec) # minute
hr  = unit("hr", 60 * Min) # hour

ft  = unit("ft") # foot
FL  = unit("FL", 100 * ft) # Flight Level (100 ft)
kft = unit("kft", 1000 * ft) # kilo-feet (nonstandard but convenient)
nmi = unit("nmi", 6076.11549 * ft) # nautical mile
km  = unit("km", nmi/1.852) # kilometer

lbm = unit("lbm") # pound mass
kg  = kilogram = unit("kg", lbm/0.45359237)  # mass

kn  = unit("kn", nmi/hr) # knots
fpm = unit("fpm", ft/Min) # ft/min (for altitude rate)

gacc = unit("gacc", 32.2 * ft/sec**2) # gravitational acceleration at sea level
lbf = unit("lbf", gacc * lbm) # pound force

unit_("nmi/min", nmi/Min)

unit_type(rad, "dimensionless (angle)")
unit_type(sec, "time")
unit_type(ft, "length")
unit_type(lbm, "mass")
unit_type(kn, "speed")
unit_type(gacc, "acceleration")
unit_type(lbf, "force")

def stop(): return raw_input() # for debugging convenience (easier to type)
