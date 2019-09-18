from __future__ import division # let 5/2 = 2.5 rather than 2
from ACPropeller import ACPropeller
from ACEngine import ACEngine
from ACPropulsion import ACPropulsion
import pylab as pyl
import numpy as npy
from scalar.units import IN, LBF, PSFC, SEC, ARCDEG, FT, OZF, RPM, HP

# Set Propeller properties
Prop = ACPropeller()
Prop.name       = 'Prop 14x5'
Prop.D          = 14.2*IN
Prop.LenDi      = [0*IN, Prop.D]
Prop.PitchAngle = 12*ARCDEG
Prop.dAlpha     = 0*ARCDEG
Prop.Solidity   = 0.0136
Prop.RD         = 3/8
Prop.AlphaStall = 14*ARCDEG
Prop.Weight     = 0.3*LBF


# Set Engine properties
Engine  = ACEngine()
Engine.name         = 'Engine'
Engine.Rbs          = 1.1
Engine.Rla          = 3.5
Engine.NumCyl       = 1
Engine.NumRev       = 1
Engine.CompRatio    = 9
Engine.Vd           = 0.607*IN**3
Engine.PistonSpeedR = 38.27*FT/SEC
Engine.BMEPR        = 50.038*LBF/IN**2
Engine.Rnmt         = 0.01
Engine.Rtmt         = 1.5
Engine.MEPtlmt      = 10.1526*LBF/IN**2
Engine.SFCmt        = 1*PSFC
Engine.Rnsfc        = 0.8
Engine.A_F          = 16
Engine.PS           = 1
Engine.LWH          = [3*IN, 2*IN, 2*IN]
Engine.Weight       = 1.5*LBF

# Set Propulsion properties
Propulsion = ACPropulsion(Prop,Engine)
Propulsion.Alt  = 0*FT
Propulsion.Vmax = 100*FT/SEC
Propulsion.nV   = 20

if __name__ == '__main__':
    
    Propulsion.Refresh()
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
    pyl.show()
    