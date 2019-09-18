#!/usr/bin/env python

"""
Copyright (c) 2006, Russell A. Paielli <http://RussP.us>
All rights reserved.

See license file for terms of license.

This file implements the SI (System Internationale) metric system of
units and many other commonly used non-metric units. A user guide is
available at <http://RussP.us/scalar.htm>. Send comments, corrections or
suggestions to Russ Paielli at Russ.Paielli@gmail.com or use the contact
link at <http://RussP.us>.

The units defined in this file are based on the inAsUnition from the
following websites:

http://physics.nist.gov/cuu/Units/units.html
http://www.bipm.org/en/si/
http://www.ex.ac.uk/cimt/dictunit/dictunit.htm
http://ts.nist.gov/ts/htdocs/230/235/appxc/appxc.htm
http://en.wikipedia.org/wiki/Category:Units_of_measure
"""

from __future__ import division # let 5/2 = 2.5 rather than 2

from os import environ as _environ

if "scalar_off" in _environ: from scalar_off import *
else: from scalar import *

# SI base units:

SEC = second     = unit("s")   # time
M = mtr = meter  = unit("m")   # length
KG = kilogram    = unit("kg")  # mass
A = amp = ampere = unit("A")   # electric current
K = kel = kelvin = unit("K")   # temperature
mol = mole       = unit("mol") # amount of substance
cd = candela     = unit("cd")  # luminous intensity
RAD = radian     = unit()      # plane angle

# Common scaled variations of base units:

ms = millisecond = unit("ms", SEC/1000)  # time
us = microsecond = unit("us", ms/1000) # time
MM = millimeter  = unit("mm", M/1000)  # length
um = micrometer  = unit("um", MM/1000) # length
cm = centimeter  = unit("cm", M/100)   # length
km = kilometer   = unit("km", 1000*M)  # length
g  = GRAM        = unit("g" , KG/1000) # mass
mg = milligram   = unit("mg", g/1000)  # mass
ug = microgram   = unit("ug", mg/1000) # mass
mA = milliampere = unit("mA", A/1000)  # electric current
uA = microampere = unit("uA", mA/1000) # electric current

# Derived units with special names and symbols:


sr  = steradian = unit("sr", RAD) # solid angle
ARCDEG = degree = unit("deg", RAD * pi/180) # plane angle
HZ  = hertz = unit("Hz", 1./SEC) # frequency
N   = newton = unit("N", KG*M/SEC**2) # force
Pa  = pascal = unit("Pa", N/M**2) # pressure
J   = joule = unit("J", N*M) # energy, work, or quantity of heat
W   = WATT = unit("W", J/SEC) # power
C   = coulomb = unit("C", A*SEC) # electric charge
V   = volt = unit("V", W/A) # electric potential
F   = farad = unit("F", C/V) # capacitance
OHM = unit("ohm", V/A) # electric resistance
S   = siemens = unit("S", A/V) # electric conductance
Wb  = weber = unit("Wb", V*SEC) # magnetic flux
T   = tesla = unit("T", Wb/M**2) # magnetic flux density
H   = henry = unit("H", Wb/A) # inductance
lm  = lumen = unit("lm", cd*sr) # luminous flux
lx  = lux = unit("lx", lm/M**2) # illuminance
Bq  = becquerel = unit("Bq", 1./SEC) # activity (of a radionuclide)
Gy  = gray = unit("Gy", J/KG) # absorbed dose, specific energy (imparted), kerma
Sv  = sievert = unit("Sv", J/KG) # dose equivalent
kat = katal = unit("kat", mol/SEC) # catalytic activity

# Common scaled variations of derived units:

kHz = kilohertz = unit("kHz", 1000 * HZ) # frequency
MHz = megahertz = unit("MHz", 1000 * kHz) # frequency
GHz = gigahertz = unit("GHz", 1000 * MHz) # frequency
kN  = kilonewton = unit("kN", 1000 * N) # force
MPa = megapascal = unit("MPa", 1000000 * Pa) # pressure
kPa = kilopascal = unit("kPa", 1000 * Pa) # pressure
hPa = hectopascal = unit("hPa", 100 * Pa) # pressure
kJ  = kilojoule = unit("kJ", 1000 * J) # energy, work, or quantity of heat
MJ  = megajoule = unit("MJ", 1000 * kJ) # energy, work, or quantity of heat
GJ  = gigajoule = unit("GJ", 1000 * MJ) # energy, work, or quantity of heat
kW  = kilowatt = unit("kW", 1000 * W) # power
MW  = megawatt = unit("MW", 1000 * kW) # power
GW  = gigawatt = unit("GW", 1000 * MW) # power
kV  = kilovolt = unit("kV", 1000 * V) # electric potential
mS  = millisiemens = unit("mS", S/1000) # electric conductance
mW  = milliwatt = unit("mW", W/1000) # power
mH  = millihenry = unit("mH", H/1000) # inductance
uH  = microhenry = unit("uH", mH/1000) # inductance
mF  = millifarad = unit("mF", F/1000) # capacitance
uF  = microfarad = unit("uF", mF/1000) # capacitance

# Output conversion for derived units expressed in terms of base units:

unit_("m**2", M**2) # square meter: area
unit_("m**3", M**3) # cubic meter: volume
unit_("m/s", M/SEC) # meter/second: speed, velocity
unit_("m/s**2", M/SEC**2) # meter/second^2: acceleration
unit_("/m", 1/M) # reciprocal meter: wave number
unit_("kg/m**3", KG/M**3) # mass density
unit_("m**3/kg", M**3/KG) # specific volume
unit_("A/m**2", A/M**2) # current density
unit_("A/m", A/M) # magnetic field strength
unit_("mol/m**3", mol/M**3) # substance concentration
unit_("cd/m**2", cd/M**2) # luminance

unit_("Pa*s", Pa*SEC) # dynamic viscosity
unit_("N*m", N*M) # moment of force
unit_("N/m", N/M) # surface tension
unit_("rad/s", RAD/SEC) # angular velocity
unit_("rad/s**2", RAD/SEC**2) # angular acceleration
unit_("W/m**2", W/M**2) # heat flux density, irradiance
unit_("J/K", J/K) # heat capacity, entropy
unit_("J/(kg*K)", J/(KG*K)) # spec. heat capacity or entropy
unit_("J/kg", J/KG) # specific energy
unit_("W/(m*K)", W/(M*K)) # thermal conductivity
unit_("J/m**3", J/M**3) # energy density
unit_("V/m", V/M) # electric field strength
unit_("C/m**3", C/M**3) # electric charge density
unit_("C/m**2", C/M**2) # electric flux density
unit_("F/m", F/M) # permittivity
unit_("H/m", H/M) # permeability
unit_("J/mol", J/mol) # molar energy
unit_("J/(mol*K)", J/(mol*K)) # molar entropy
unit_("C/kg", C/KG) # exposure (to x and gamma rays)
unit_("Gy/s", Gy/SEC) # absorbed dose rate
unit_("W/sr", W/sr) # radiant intensity
unit_("W/(m**2*sr)", W/sr) # radiance
unit_("kat/m**3", kat/M**3) # catalytic (activity) concentration

unit_("deg/s", ARCDEG/SEC ) # angular velocity
unit_("rad/s", RAD/SEC    ) # angular velocity

# other common units (including non-metric):

MIN = minute = unit("min", 60 * SEC) # time
HR  = hour = unit("hr", 60 * MIN) # time

FT  = feet = unit("ft", 0.3048 * M) # length
FL  = FlightLevel = unit("FL", 100 * FT) # length (for aircraft altitude)
kft = kilofeet = unit("kft", 1000 * FT) # length (for aircraft altitude)
IN  = inch = unit("in", FT/12.) # length
nmi = nmile = unit("nmi", 1852 * M) # length
mi  = mile = unit("mi", 1609.344 * M) # US/UK statute mile: length

unit_("in**2", IN**2)
unit_("in**3", IN**3)
unit_("ft**2", FT**2)
unit_("ft**3", FT**3)
unit_("mi**2", mi**2)
unit_("mi**3", mi**3)
unit_("nmi**2", nmi**2)
unit_("nmi**3", nmi**3)

LBM = unit( "lbm", 0.45359237 * KG) # pound mass
OZM = unit( "ozm", LBM / 16       ) # ounce mass
LBF = unit( "lbf", 4.448222 * N   ) # pound force
OZF = unit( "ozf", LBF / 16       ) # ounce force
dyn = dyne = unit("dyn", 1e-5 * N ) # force
SLUG = unit("slug", 14.59390 * KG ) # mass

kn  = knot = unit("kn", nmi/HR) # knots: speed
kph = unit("kph", km/HR) # kilometers per hour: speed
mph = unit("mph", mi/HR) # miles per hour: speed
fpm = unit("fpm", FT/MIN) # speed (for altitude rate)
fps = unit("fps", FT/SEC) # speed

gacc = unit("gacc", 9.80665 * M/SEC**2) # gravitational accel at sea level

BAR = unit("bar", 100 * kPa) # bar: pressure
ATM = unit("atm", 101.325 * kPa) # atmospheres: pressure
PSI = unit("psi", LBF/IN**2) # pounds per square inch: pressure
inHg = unit("inHg", 3386.388 * Pa) # inches of Mercury: pressure

HP = unit("hp", 33000.*FT*LBF/MIN) # horsepower

Wh  = unit("Wh", W * HR) # Watt-hour: energy
kWh = unit("kWh", 1000 * Wh) # kilowatt-hour: energy
Btu = unit("Btu", 1055.056 * J) # British thermal unit: energy (Intl Table)
erg = unit("erg", 1e-7 * J) # erg: energy
eV  = unit("eV", 1.602176e-19 * J) # electronvolt: energy
cal = unit("cal", 4.185 * J) # gram ("small") calorie: energy

RPM = unit("rpm", 2*pi*RAD/MIN) # revolutions per minute: angular velocity
cfm = unit("cfm", FT**3/MIN) # cubic feet per minute: flow rate

acre = unit("acre", 43560 * FT**2) # acre: area
ha = hectare = unit("ha", 1e4 * M**2) # area

GAL = gallon = unit("gal", 231*IN**3) # US gallon: volume (usually of liquid)
L = liter = unit("L", 1000 * cm**3) # volume (usually of liquid)

degF = unit("degF", K / 1.8) # (delta) degrees Fahrenheit: temperature
degR = unit("degR", K / 1.8) # (delta) degrees Rankine: temperature
degC = unit("degC", K) # (delta) degrees Celsius: temperature

Oe = oersted = unit("Oe", 79.577 * A/M) # magnetic field strength
fc = footcandle = unit("fc", lm/FT**2) # illuminance (light intensity)

Ah = amphour = unit("Ah", A*HR)
mAh = milliamphour = unit("mAh", A*HR/1000)


# Fuel Consumption
TSFC = unit( "tsfc", LBM/(HR*LBF) ) # Thrust specific fuel consumption
PSFC = unit( "psfc", LBM/(HR*HP) )  # Power specific fuel consumption


unit_("nmi/min", nmi/MIN)
unit_("ft*lbf", FT*LBF)
unit_("1/rad", 1./RAD)
unit_("ft/min", FT/MIN)
unit_("lbf/in", LBF/IN)
unit_("lbf/in**2", LBF/IN**2)
unit_("lbf/ft", LBF/FT)
unit_("lbf/ft**2", LBF/FT**2)

unit_("ft/s", FT/SEC) # meter/second: speed, velocity
unit_("slug*ft**2", SLUG*FT**2) # moment of inertia
unit_("in*ozf", IN*OZF)
unit_("in*lbf", IN*LBF)
unit_("lbf*in", LBF*IN)
unit_("g/in**3", GRAM/IN**3)


unit_type(SEC, "time")
unit_type(M, "length")
unit_type(KG, "mass")
unit_type(A, "electric current")
unit_type(K, "temperature")
unit_type(mol, "amount of substance")
unit_type(cd, "luminous intensity")

unit_type(M**2, "area")
unit_type(M**3, "volume")
unit_type(M/SEC, "speed")
unit_type(M/SEC**2, "acceleration")
unit_type(1/M, "wave number")
unit_type(KG/M**3, "mass density")
unit_type(M**3/KG, "specific volume")
unit_type(A/M**2, "current density")
unit_type(A/M, "magnetic field strength")
unit_type(mol/M**3, "substance concentration")
unit_type(cd/M**2, "luminance")

unit_type(RAD, "dimensionless (or angle)")
unit_type(HZ, "frequency")
unit_type(N, "force")
unit_type(Pa, "pressure")
unit_type(J, "energy (or moment of force)")
unit_type(W, "power")
unit_type(C, "electric charge")
unit_type(V, "electric potential")
unit_type(F, "capacitance")
unit_type(OHM, "electric resistance")
unit_type(S, "electric conductance")
unit_type(Wb, "magnetic flux")
unit_type(T, "magnetic flux density")
unit_type(H, "inductance")
unit_type(lm, "luminous flux")
unit_type(lx, "illuminance")
unit_type(Bq, "activity (of a radionuclide)")
unit_type(Sv, "dose equivalent")
unit_type(kat, "catalytic activity")

unit_type(Pa*SEC, "dynamic viscosity")
unit_type(N/M, "surface tension")
unit_type(RAD/SEC, "angular speed")
unit_type(RAD/SEC**2, "angular acceleration")
unit_type(W/M**2, "heat flux density")
unit_type(J/K, "heat capacity or entropy")
unit_type(J/(KG*K), "specific heat capacity or entropy")
unit_type(J/KG, "specific energy")
unit_type(W/(M*K), "thermal conductivity")
unit_type(J/M**3, "energy density")
unit_type(V/M, "electric field strength")
unit_type(C/M**3, "electric charge density")
unit_type(C/M**2, "electric flux density")
unit_type(F/M, "permittivity")
unit_type(H/M, "permeability")
unit_type(J/mol, "molar energy")
unit_type(J/(mol*K), "molar entropy")
unit_type(C/KG, "exposure (x and gamma rays)")
unit_type(Gy/SEC, "absorbed dose rate")
unit_type(kat/M**3, "catalytic (activity) concentration")
