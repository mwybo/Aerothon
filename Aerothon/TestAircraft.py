from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import M, FT, IN, ARCDEG, RAD, LBF, SEC, KG
from ACWingWeight import ACSolidWing
from DefaultMaterialsLibrary import PinkFoam
from ACAircraft import ACTailAircraft
from TestFuselage import Fuselage
from TestPropulsion import Propulsion
import pylab as pyl

aircraft = ACTailAircraft()

aircraft.SetFuselage(Fuselage)
aircraft.SetPropulsion(Propulsion)

aircraft.TotalWeight = 25*LBF
aircraft.EmptyWeight = 5.256*LBF

aircraft.TippingAngle  = 10*ARCDEG
aircraft.RotationAngle = 10*ARCDEG

aircraft.Alpha_Groundroll = 0*ARCDEG
    
aircraft.CMSlopeAt   = (1.5 * ARCDEG, 2.5 * ARCDEG)
aircraft.CLSlopeAt   = (6 * ARCDEG, 7 * ARCDEG)
aircraft.CLHTSlopeAt = (1.5 * ARCDEG, 2.5 * ARCDEG) 
aircraft.DWSlopeAt   = (1 * ARCDEG, 2 * ARCDEG)

aircraft.StaticMargin   = 0.1
aircraft.Alpha_Zero_CM  = 3 *ARCDEG
aircraft.Xnp_Init       = 40.5 *IN

aircraft.Wing.Lift_LO       = 26.799 * LBF
aircraft.Wing.V_Stall       = 39.2 * FT/SEC
aircraft.Wing.Alt_LO        = 600 * FT
aircraft.Wing.ClSlopeAt     = (6 * ARCDEG, 7 * ARCDEG)

#aircraft.Wing.S             = 900 * IN**2

#aircraft.Wing.AR            = 16
aircraft.Wing.b             = 10. * FT
aircraft.Wing.TR            = [1, 0.7, 0.7]
aircraft.Wing.Fb            = [0.3, 0.8, 1]
#aircraft.Wing.Gam           = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
aircraft.Wing.CEdge         = 'LE'
aircraft.Wing.ConstUpper    = True
aircraft.Wing.Airfoil       = 'S1223'
aircraft.Wing.FullWing      = True
aircraft.Wing.LLFile        = 'LLTPolar.txt'

#
# Make the wing a solid foam wing
#
aircraft.Wing.SetWeightCalc(ACSolidWing)
aircraft.Wing.WingWeight.AddSpar("MainSpar",1*IN,0.5*IN)
aircraft.Wing.WingWeight.MainSpar.SparMat           = PinkFoam.copy()
aircraft.Wing.WingWeight.SkinMat.AreaForceDensity   = 0.001*LBF/IN**2
aircraft.Wing.WingWeight.WingMat                    = PinkFoam.copy()

aircraft.Wing.Aileron.Fc = 0.25
aircraft.Wing.Aileron.Fb = 0.42
aircraft.Wing.Aileron.Ft = 0.2

aircraft.HTail.Airfoil  = 'NACA0012'
aircraft.HTail.AR       = 3
aircraft.HTail.TR       = 0.5
aircraft.HTail.S        = 100 * IN**2
#aircraft.HTail.L        = 30 * IN
aircraft.HTail.o_eff    = 0.96
aircraft.HTail.VC       = 0.7
aircraft.HTail.FullWing = True
aircraft.HTail.DWF      = 2.0 #Main wing Down wash factor

aircraft.HTail.Controls.Elevator.Fc = 0.3

#Set the sweep about the elevator hinge
aircraft.HTail.SweepFc  = 1.0 - aircraft.HTail.Controls.Elevator.Fc


aircraft.VTail.Airfoil      = 'NACA0012'
aircraft.VTail.VC           = 0.06
aircraft.VTail.AR           = 1.6
aircraft.VTail.TR           = 0.7
aircraft.VTail.Axis         = (0, 1)
#aircraft.VTail.L            = 51.572 * IN
aircraft.VTail.S            = 69 * IN**2

aircraft.VTail.Controls.Rudder.Fc = 0.4

#Set the sweep about the rudder hinge
aircraft.VTail.SweepFc      = 1.0 - aircraft.VTail.Controls.Rudder.Fc


# Aircraft Landing Gear
MainGear = aircraft.MainGear
MainGear.Theta = 45*ARCDEG
MainGear.GearHeight, MainGear.StrutW, MainGear.StrutH = 5 * IN, 0.2 * IN, 0.1 * IN
MainGear.X = [5 * IN, 2 *IN, MainGear.GearHeight]

NoseGear = aircraft.NoseGear
NoseGear.Theta = 0*ARCDEG
NoseGear.GearHeight, NoseGear.StrutW, NoseGear.StrutH = 5 * IN, 0.2 * IN, 0.1 * IN
NoseGear.X = [0 * IN, 0 * IN, MainGear.GearHeight]


### ACCT = ACControls(aircraft)

if __name__ == '__main__':
    
    #Xcg = aircraft.Wing.Xac() + 2. * IN
    #del_e = 0 * ARCDEG
    #alpha2dw = 5 * ARCDEG
    #print 'CM', aircraft.CM(alpha2dw, del_e)
    #print 'dCM_da', aircraft.dCM_da(del_e)
    #print 'Wing AC', aircraft.Wing.Xac()
    #print 'Neutral point', aircraft.Xnp()
    #print 'Center of Gravity', aircraft.Xcg()
    #print 'HTail i', aircraft.HTail.i
    print 'Aircraft Xnp   :', aircraft.Xnp()
    print 'Aircraft Xcg   :', aircraft.Xcg()
    print 'Aircraft CM    :', aircraft.CM(15*ARCDEG, del_e = 10*ARCDEG)
    print 'Aircraft dCM_da:', aircraft.dCM_da(0*ARCDEG, aircraft.Xnp())
    print
    print 'Wing Area      :', aircraft.Wing.S
    print 'Wing MAC       :', aircraft.Wing.MAC()
    print 'Wing dCl_da    :', aircraft.Wing.dCL_da()
    print 'Wing Xac       :', aircraft.Wing.Xac()
    print 'Wing Weight    :', aircraft.Wing.Weight
    print
    print 'Horiz Area     :', aircraft.HTail.S
    print 'Horiz Length   :', aircraft.HTail.L
    print 'Horiz iht      :', aircraft.HTail.i
    print 'Horiz dCl_da   :', aircraft.HTail.dCL_da()
    print
    print 'Vert Area      :', aircraft.VTail.S
    print 'Vert Span      :', aircraft.VTail.b
    
    #alphas = aircraft.Wing.AlphaRange()
    #del_e_trim = aircraft.del_e_trim(alphas)
    
    #pyl.figure(3)
    #pyl.plot(alphas / (ARCDEG),del_e_trim / (ARCDEG))
    
    #print
    #print 'Aircraft del_e_trim:', aircraft.del_e_trim(5*ARCDEG)
    
    #aircraft.WriteAVLMainWing('AVLMainWing.avl')
    #subp.call(['avl26.exe', 'TestAVL.avl < AVLCommands.dat'])
    
    #aircraft.Wing.LLTPlots(3)
    #aircraft.Wing.Draw2DPolar(4)
    #aircraft.HTail.Draw2DPolar(3)
    #aircraft.DrawCMPlot(2, (-10*ARCDEG, -5*ARCDEG, 0*ARCDEG, +5*ARCDEG, +10 * ARCDEG), (+0.5 * IN, -0.5 * IN))
    aircraft.Draw()
    pyl.show()
    
