
from scalar.units import IN, LBF
from ACFuselage import ACFuselage

Fuselage = ACFuselage()
Fuselage.name = "Fuselage"

Fuselage.AddSection('Nose'     , 7*IN,  1)
Fuselage.AddSection('PyldBay'  , 10*IN, 1)
Fuselage.AddSection('TailTaper', 15*IN)
Fuselage.AddSection('Tail')

#Fuselage.Nose.Length = 10*IN
Fuselage.Nose.FrontBulk.Width = 3*IN
Fuselage.Nose.FrontBulk.Height = 3*IN
Fuselage.Nose.FrontBulk.Weight = 0.5*LBF
Fuselage.Nose.Weight = 0.2*LBF


Fuselage.PyldBay.FrontBulk.Width = 5*IN
Fuselage.PyldBay.FrontBulk.Height = 5*IN
Fuselage.PyldBay.SkinMat.AreaForceDensity = 0.01*LBF/IN**2

Fuselage.PyldBay.BackBulk.Width = 5*IN
Fuselage.PyldBay.BackBulk.Height = 5*IN

Fuselage.TailTaper.Align = 1

Fuselage.Tail.BackBulk.X[0] = 50*IN
Fuselage.Tail.BackBulk.X[2] = 0*IN
Fuselage.Tail.Weight = 0.7*LBF

Fuselage.Nose.AddComponent    (     "Battery"   , 0.5*LBF , (0.25*IN,1.5*IN,1*IN)    , "Right"  , (0.2,0.5,0.5)   )
#Fuselage.PyldBay.AddComponent (     "Payload"   , 15.*LBF , (5*IN,3.5*IN,3*IN)       ,  None    , (0.2, 0.5, 0.2)   )
Fuselage.Nose.AddComponent    (     "FuelTank"  , 0.5*LBF , (2.5*IN,2*IN,1.25*IN)    , "Back"   , (0.75, 0.5, 0.7) )
Fuselage.Nose.AddComponent    ("NoseWheelServo" , 0.1*LBF , (.5*IN,1*IN,1*IN)        , "Bottom" , (0.6, 0.2, 0) )

Fuselage.XcgSection = Fuselage.PyldBay
Fuselage.TailBulk = Fuselage.Tail.BackBulk

