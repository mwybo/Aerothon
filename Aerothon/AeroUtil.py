"""
A collection of useful atmospheric and otherwise related functions
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import FT, degR, ATM, LBM, LBF, KG, M, SEC
import numpy as npy
from pylab import plot, show, figure, xlabel, legend
import cmath as math

#
# Standard Atmosphere Functions
#
BTU = 25036.856 * LBM*FT**2/SEC**2

# Sea level conditions
Tsl   = 518.69 * degR               # Temperature
Psl   = 1 * ATM                     # Pressure
rhosl = 0.076 * LBM / FT**3         # Density
musl  = 1.783 * 10**-5 * KG/(M*SEC) # Viscosity
gam   = 1.4                         # Ratio of specific heats
Cp    = 0.238 *BTU/(LBM*degR)       # Specific heat at constant pressure
R     = (gam-1)/gam*Cp              # Gas constant

# Temperature
def TempAlt(h):
    """
    Computes temperature

    Inputs:
        h - Altitude
    """
    ah = h / FT
    Alt  = npy.array([0     , 36152, 82346, 155348, 175344, 250000, 315000])
    Temp = npy.array([518.69, 390.0, 390.0,  508.8,  508.8, 349.67, 349.67])

    return npy.interp(ah, Alt, Temp) * degR

# Pressure
def PressAlt(h):
    """
    Computes pressure

    Inputs:
        h - Altitude
    """
    ah = h / FT
    if ah >= -15000 and ah < 36089:
        P = abs((1 - ah/(145442.3)) ** 5.25591)
    elif ah >= 36089 and ah < 65617:
        P = abs(0.22336 * math.exp((36089.24 - ah)/20805.722))
    elif ah >= 65617 and ah < 104987:
        P = abs(0.0540322 * ((645175.9 + ah)/710792.7) ** -34.1634)
    elif ah >= 104987:
        P = abs(0.00856649 * ((162928.62 + ah)/267915.5) ** -12.2012)

    return P * ATM

# Density
def DensityRatioAlt(h):
    """
    Computes density at altitude

    Inputs:
        h - Altitude
    """
    T = TempAlt(h)/Tsl
    P = PressAlt(h)/Psl

    return P/T

def DensAlt(h):
    """
    Computes density at altitude

    Inputs:
        h - Altitude
    """
    return DensityRatioAlt(h)*rhosl

def STDCorrection(P, T):
    """
    Computes the correction for standard day

    Inputs:
        P - Pressure
        T - Temperature
    """
    return (Psl/Tsl)/(P/T)

# Molecular viscosity
def Mu(T):
    """
    Computes molecular viscosity based on Sutherlands Law
    This version has a 2% error between 300 - 3420 degR

    Inputs:
        T - Temperature
    """
    mu0 = 2.27*10**-8 * LBF*SEC/(FT**2 * degR**(1./2.))
    return mu0 * (T**(3./2.)/(T + 198.6*degR))

# Reynolds number
def Re(h, v, l):
    """
    Computes Reynolds number

    Inputs:
        h - Altitude
        v - Velocity
        l - Reference Length
    """
    return DensAlt(h) * v * l/Mu(TempAlt(h))

# Dynamic pressure
def q(h, v):
    """
    Computes dynamic pressure

    Inputs:
        h - Altitude
        v - Velocity
    """
    return DensAlt(h) * v **2 / 2

# Speed of sound
def SpeedOfSound(h):
    """
    Computes speed of sound at a given altitude

    Inputs:
        h - Altitude
    """
    T = TempAlt(h)
    return (gam*R*T)**0.5

# Mach number
def Mach(h, v):
    """
    Computes speed of sound at a given altitude

    Inputs:
        h - Altitude
        v - Velocity
    """
    a = SpeedOfSound(h)
    return v/a

if __name__ == '__main__':
    x = []
    T = []
    P = []
    rho = []
    for i in range(0, 315000, 1000):
        x.append(i)
        T.append((TempAlt(i*FT)/Tsl))
        P.append((PressAlt(i*FT)/Psl))
        rho.append((DensAlt(i*FT)/rhosl))

    Alt = 600 * FT
    print 'Quantities at ', Alt
    print 'Temperature:', TempAlt(Alt)
    print 'Pressure:', PressAlt(Alt)
    print 'Density:', DensAlt(Alt)
    print 'Dynamic viscosity:', Mu(TempAlt(600 * FT))
    print 'Reynolds Number:', Re(Alt, 1 * FT/SEC, 1 * FT)
    
    Alt = 0 * M
    print 'Temperature:', TempAlt(Alt)
    print 'Reynolds Number:', Re(Alt, 90 * M/SEC, 0.07 * M)
    
    Alt = 5000 * M
    print 'Temperature:', TempAlt(Alt)
    print 'Reynolds Number:', Re(Alt, 90 * M/SEC, 0.07 * M)
    
    figure(1)
    plot(x,T,x,P,x,rho)
    xlabel("Altitude (ft)")
    legend(["T","P","rho"])
    show()    

