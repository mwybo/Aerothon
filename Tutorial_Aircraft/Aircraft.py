from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import M, FT, IN, ARCDEG, RAD, LBF, SEC, KG, SLUG
from scalar.units import AsUnit
from Aerothon.ACAircraft import ACTailAircraft
from Aerothon.ACWingWeight import ACSolidWing
from Aerothon.DefaultMaterialsLibrary import Monokote, PinkFoam, Balsa
from Fuselage import Fuselage
from Propulsion import Propulsion
from Wing import Wing
import pylab as pyl

#
# Create the Aircraft
#
Aircraft = ACTailAircraft()
Aircraft.name = 'Tutorial Aircraft'

# 
# Assign the already generated parts
#
Aircraft.SetFuselage(Fuselage)
Aircraft.SetPropulsion(Propulsion)
Aircraft.SetWing(Wing)

#
# Aircraft Properties
#
Aircraft.TotalWeight = 28*LBF

Aircraft.TippingAngle     = 10*ARCDEG
Aircraft.RotationAngle    = 10*ARCDEG
Aircraft.Alpha_Groundroll = 0*ARCDEG # incidence on wing
    
Aircraft.CMSlopeAt   = (0 * ARCDEG, 3 * ARCDEG) #plotting
Aircraft.CLSlopeAt   = (6 * ARCDEG, 7 * ARCDEG) #plotting
Aircraft.CLHTSlopeAt = (0 * ARCDEG, 9 * ARCDEG) #plotting
Aircraft.DWSlopeAt   = (7 * ARCDEG, 8 * ARCDEG)#plotting

Aircraft.Alpha_Zero_CM  = 0 * ARCDEG
Aircraft.StaticMargin   = 0.1

#
# Maximum velocity for plotting purposes
#
Aircraft.VmaxPlt = 100*FT/SEC

###############################################################################
#
# Tail surfaces
#
###############################################################################

#==============================================================================
# Horizontal tail
#
HTail = Aircraft.HTail
HTail.Airfoil  = 'NACA0012'
HTail.AR       = 3.5
HTail.TR       = 0.5
HTail.S        = 200 * IN**2
#HTail.L        = 30 * IN
HTail.VC       = 0.65
HTail.FullWing = True
HTail.DWF      = 2.0 #Main wing Down wash factor (Between 1.0 and 2.0)
HTail.Inverted = False

#
# Elevator properties
#
HTail.Elevator.Fc = 0.35
HTail.Elevator.Fb = 1.0
HTail.Elevator.Ft = 0.0

HTail.Elevator.Servo.Fc  = 0.3
HTail.Elevator.Servo.Fbc = 0.1

#Set the sweep about the elevator hinge
HTail.SweepFc  = 1.0 - HTail.Elevator.Fc

#
# Structural properties
#
HTail.SetWeightCalc(ACSolidWing)
HTail.WingWeight.AddSpar("HMainSpar",0.25*IN,0.5*IN)
HTail.WingWeight.HMainSpar.SparMat = Balsa.copy()
HTail.WingWeight.SkinMat                    = Monokote.copy()
HTail.WingWeight.WingMat                    = PinkFoam.copy()


#==============================================================================
# Vertical tail
#
VTail = Aircraft.VTail
VTail.Airfoil = 'NACA0012'
VTail.VC      = 0.05
VTail.AR      = 1.6
VTail.TR      = 0.7
VTail.Axis    = (0, 1)
#VTail.L       = 51.572 * IN
VTail.S       = 69 * IN**2
#VTail.b       = 10 * IN

#
# Rudder properties
#
VTail.Rudder.Fc = 0.4
VTail.Rudder.Servo.Fc  = 0.3
VTail.Rudder.Servo.Fbc = 0.1

#Set the sweep about the rudder hinge
VTail.SweepFc = 1.0 - VTail.Rudder.Fc

#
# Structural properties
#
VTail.SetWeightCalc(ACSolidWing)
VTail.WingWeight.AddSpar("VMainSpar",0.25*IN,0.25*IN)
VTail.WingWeight.VMainSpar.SparMat = Balsa.copy()
VTail.WingWeight.SkinMat                    = Monokote.copy()
VTail.WingWeight.WingMat                    = PinkFoam.copy()


###############################################################################
#
# Landing Gear
#
###############################################################################

MainGear = Aircraft.MainGear
MainGear.Theta        = 30*ARCDEG
MainGear.GearHeight   = 6   * IN
MainGear.StrutW       = 0.2 * IN
MainGear.StrutH       = 0.1 * IN
MainGear.WheelDiam    = 2.5 * IN
MainGear.X[1]         = 2.5 * IN
MainGear.Strut.Weight = 0.2*LBF
MainGear.Wheel.Weight = 0.1*LBF

NoseGear = Aircraft.NoseGear
NoseGear.StrutW    = 0.1 * IN
NoseGear.StrutH    = 0.1 * IN
NoseGear.WheelDiam = 2.5 * IN
NoseGear.Strut.Weight = 0.1*LBF
NoseGear.Wheel.Weight = 0.1*LBF


if __name__ == '__main__':
    
    TCG = Aircraft.CG()
    ACG = Aircraft.Fuselage.AircraftCG()
    dCG = TCG - ACG
    
    print 'Aircraft Xnp   :', AsUnit( Aircraft.Xnp(), 'in' )
    print 'Aircraft Xcg   :', AsUnit( Aircraft.Xcg(), 'in' )
    print 'Aircraft CM    :', Aircraft.CM(15*ARCDEG, del_e = 10*ARCDEG)
    print 'Aircraft dCM_da:', AsUnit( Aircraft.dCM_da(), '1/rad' )
    print
    print 'Wing Area      :', AsUnit( Aircraft.Wing.S, 'in**2' )
    print 'Wing MAC       :', AsUnit( Aircraft.Wing.MAC(), 'in' )
    print 'Wing dCl_da    :', AsUnit( Aircraft.Wing.dCL_da(), '1/rad' )
    print 'Wing Xac       :', AsUnit( Aircraft.Wing.Xac(), 'in' )
    print
    print 'Horiz Area     :', AsUnit( Aircraft.HTail.S, 'in**2' )
    print 'Horiz Length   :', AsUnit( Aircraft.HTail.L, 'in' )
    print 'Horiz iht      :', AsUnit( Aircraft.HTail.i, 'deg' )
    print 'Horiz dCl_da   :', AsUnit( Aircraft.HTail.dCL_da(), '1/rad' )
    print 'Horiz Design CL:', Aircraft.GetHT_Design_CL()
    print
    print 'Wing  Weight   :', AsUnit( Aircraft.Wing.Weight, 'lbf' )
    print 'Horiz Weight   :', AsUnit( Aircraft.HTail.Weight, 'lbf' )
    print 'Vert  Weight   :', AsUnit( Aircraft.VTail.Weight, 'lbf' )
    print
    print 'Vert Area      :', AsUnit( Aircraft.VTail.S, 'in**2' )
    print 'Vert Span      :', AsUnit( Aircraft.VTail.b, 'in' )
    print
    print 'True CG        :', AsUnit( TCG[0], 'in' )
    print 'Desired CG     :', AsUnit( ACG[0], 'in' )
    print 'delta CG       :', AsUnit( dCG[0], 'in' )
    print 'Empty   Weight :', AsUnit( Aircraft.EmptyWeight, 'lbf' )
    print 'Payload Weight :', AsUnit( Aircraft.PayloadWeight(), 'lbf' )
    print
    print 'Propulsion MOI :', AsUnit( Aircraft.Propulsion.MOI(), 'slug*ft**2', '%3.5f' )
    print 'Wing       MOI :', AsUnit( Aircraft.Wing.MOI(), 'slug*ft**2', '%3.5f' )
    print 'HTail      MOI :', AsUnit( Aircraft.HTail.MOI(), 'slug*ft**2', '%3.5f' )
    print 'VTail      MOI :', AsUnit( Aircraft.VTail.MOI(), 'slug*ft**2', '%3.5f' )
    print 'Aircraft   MOI :', AsUnit( Aircraft.MOI(), 'slug*ft**2', '%3.5f' )
    print
    print 'Aircraft Groundroll :', AsUnit( Aircraft.Groundroll(), 'ft' )
    
#    Aircraft.PlotCLCMComponents(fig = 7, del_es = (-10*ARCDEG, -5*ARCDEG, 0*ARCDEG, +5*ARCDEG, +10 * ARCDEG))
    Aircraft.PlotTailLoad(fig=6)
#    Aircraft.PlotCMPolars(5, (-10*ARCDEG, -5*ARCDEG, 0*ARCDEG, +5*ARCDEG, +10 * ARCDEG), (+0.5 * IN, -0.5 * IN))
    Aircraft.PlotPolarsSlopes(fig=4)
#    Aircraft.PlotDragBuildup(fig=3)
    Aircraft.PlotPropulsionPerformance(fig=2)
    Aircraft.Draw()
    
    Aircraft.WriteAVLAircraft('AVLAircraft.avl')
    pyl.show()
    
