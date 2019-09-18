from __future__ import division # let 5/2 = 2.5 rather than 2
from Aerothon.ACPropeller import ACPropeller
from Aerothon.ACEngine import ACEngine
from Aerothon.ACPropulsion import ACPropulsion
import numpy as npy
from scalar.units import IN, LBF, PSFC, SEC, ARCDEG, FT, OZF, RPM, HP
from scalar.units import AsUnit

# Set Propeller properties
Prop = ACPropeller()
Prop.D          = 13.5*IN # Diameter
Prop.Thickness  = .5*IN # Thickness at the hub... just for drawing purposes
Prop.PitchAngle = 14*ARCDEG # Pitch angle.. aka \Beta
Prop.dAlpha     = 0*ARCDEG # Difference between measured alpha and zero lift alpha
Prop.Solidity   = 0.0136 # Proportional to the blade disk area, similar to the activity factor (AreaBlades/(2*D**2))
Prop.RD         = 3/8 # The location on the profile chord where the PitchAngle is defined (default 3/8) Von Mises 306
Prop.AlphaStall = 12*ARCDEG # Stall angle of attack
Prop.Weight     = 3/32*LBF # Weight

# Use these parameters to match test data if need be.
# Prop.CLSlope    = .078/ARCDEG  #- 2D airfoil lift slope
# Prop.CDCurve    = 2.2          #- 2D curvature of the airfoil drag bucket
# Prop.CDp        = .02          #- Parasitic drag

# Set Engine properties - glow engine... see 2015 files for setting up electric motor
Engine  = ACEngine()
Engine.Rbs          = 1.1
Engine.Rla          = 3.5
Engine.NumCyl       = 1
Engine.NumRev       = 1
Engine.CompRatio    = 9
Engine.Vd           = 0.607*IN**3
Engine.PistonSpeedR = 38.27*FT/SEC
Engine.MEPtlmt      = 10.1526*LBF/IN**2
Engine.SFCmt        = 1*PSFC
Engine.A_F          = 16
Engine.PS           = 1
Engine.LWH          = [4*IN, 2.5*IN, 4*IN]
Engine.Xat          = [0, 0.5, 0.5]
Engine.Weight       = 1.5*LBF

Engine.NoseCone.LenDi = [2*IN,1.125*IN]
Engine.Muffler.LenDi = [6*IN,1.5*IN]
Engine.DrawDetail   = True

#
# Curvefitting parameters
#
Engine.BMEPR        = 50.038*LBF/IN**2
Engine.Rnmt         = 0.01
Engine.Rtmt         = 1.5
Engine.Rnsfc        = 0.8


# Set Propulsion properties
Propulsion = ACPropulsion(Prop,Engine)
Propulsion.Alt  = 0*FT
Propulsion.Vmax = 100*FT/SEC
Propulsion.nV   = 20

if __name__ == '__main__':
    import pylab as pyl
    
    #
    # This data has been corrected for standard day
    #
    #            RPM,        Torque,         Power
    TestData = [(7000  *RPM, 107.35 *IN*OZF, 0.75 *HP),
                (8500  *RPM, 104.24 *IN*OZF, 0.88 *HP),
                (9500  *RPM, 101.13 *IN*OZF, 0.95 *HP),
                (9700  *RPM, 98.02  *IN*OZF, 0.94 *HP),
                (10600 *RPM, 102.69 *IN*OZF, 1.08 *HP),
                (10800 *RPM, 98.02  *IN*OZF, 1.05 *HP),
                (11000 *RPM, 101.13 *IN*OZF, 1.1  *HP),
                (11200 *RPM, 99.57  *IN*OZF, 1.11 *HP),
                (11600 *RPM, 98.02  *IN*OZF, 1.13 *HP),
                (12900 *RPM, 93.35  *IN*OZF, 1.19 *HP),
                (13300 *RPM, 91.79  *IN*OZF, 1.21 *HP),
                (13600 *RPM, 91.79  *IN*OZF, 1.24 *HP),
                (14600 *RPM, 88.68  *IN*OZF, 1.28 *HP),
                (15600 *RPM, 79.35  *IN*OZF, 1.23 *HP),
                (16300 *RPM, 77.76  *IN*OZF, 1.26 *HP)]
    
    Engine.TestData = TestData
    Engine.PlotTestData()

    print "Static Thrust :", AsUnit( Propulsion.T(0*FT/SEC), 'lbf' )
    
    V = npy.linspace(0,100,30)*FT/SEC
    N = npy.linspace(1000,20000,30)*RPM
    #Propulsion.PlotMatched(V, fig = 2)
    V2 = npy.linspace(0,100,5)*FT/SEC
    Propulsion.PlotTPvsN(N, V2, fig=3)
    
    Propulsion.Draw(fig=4)
    
    pyl.show()
    