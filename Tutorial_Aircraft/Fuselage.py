from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import IN, LBF, SLUG, FT, OZM, OZF
from scalar.units import AsUnit
from Aerothon.ACFuselage import ACFuselage
from Aerothon.DefaultMaterialsLibrary import Monokote

Fuselage = ACFuselage()

Fuselage.AddSection('Nose'      , 6*IN  , 2)
Fuselage.AddSection('PyldBay'   , 10*IN , 1)
Fuselage.AddSection('Tail')

#
# Size the engine fire wall
#
Fuselage.Nose.FrontBulk.Width  = 2.7*IN
Fuselage.Nose.FrontBulk.Height = 2.7*IN
Fuselage.Nose.FrontBulk.Weight = 0.1*LBF
Fuselage.Nose.Weight           = 0.2*LBF

#
# Size the payload bay
#
Fuselage.PyldBay.FrontBulk.Width  = 5*IN
Fuselage.PyldBay.FrontBulk.Height = 5*IN
Fuselage.PyldBay.BackBulk.Width   = 5*IN
Fuselage.PyldBay.BackBulk.Height  = 5*IN
Fuselage.PyldBay.SkinMat.AreaForceDensity = 0.001*LBF/IN**2

#
# Give a dummy position of the last bulkhead for now
#
Fuselage.Tail.BackBulk.X[0]   = 60*IN
Fuselage.Tail.BackBulk.X[2]   = 0*IN
Fuselage.Tail.BackBulk.Width  = 0.616*IN
Fuselage.Tail.BackBulk.Height = 0.616*IN

#
# Add some components to the nose section
#
Fuselage.Nose.AddComponent    (     "Battery"   , 0.05*LBF, (0.25*IN,1.5*IN,1*IN)    , "Right"  , (0.2 , 0.5, 0.5) )
Fuselage.Nose.AddComponent    (     "FuelTank"  , 0.1*LBF , (2.5*IN,2*IN,1.25*IN)    , "Back"   , (0.75, 0.5, 0.7) )
Fuselage.Nose.AddComponent    ("NoseWheelServo" , 0.1*LBF , (.5*IN,1*IN,1*IN)        , "Bottom" , (0.6 , 0.2, 0.0) )

#
# Define which section contains the CG of the aircraft
#
Fuselage.XcgSection = Fuselage.PyldBay
Fuselage.XcgSecFrac = 0.5

#
# Determine which bulkhead should be set by the horizontal tail
#
Fuselage.TailBulk = Fuselage.Tail.BackBulk

if __name__ == '__main__':
    import pylab as pyl
    
    print 'Nose      Weight :', Fuselage.Nose.Weight
    print 'PyldBay   Weight :', Fuselage.PyldBay.Weight
    print 'Tail      Weight :', Fuselage.Tail.Weight
    
    print 'Fuselage Weight    :', Fuselage.Weight
    print 'Fuselage MOI       :', AsUnit( Fuselage.MOI(), 'slug*ft**2' )
    print 'Fuselage CG        :', AsUnit( Fuselage.CG(), 'in' )
    print 'Fuselage Desired CG:', AsUnit( Fuselage.AircraftCG(), 'in' )
    
    
    Fuselage.Draw()
    pyl.show()