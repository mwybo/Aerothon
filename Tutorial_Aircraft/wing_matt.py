from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import LBF, SEC, ARCDEG, FT, IN, SLUG
from scalar.units import AsUnit
from Aerothon.ACWing import ACMainWing
from Aerothon.DefaultMaterialsLibrary import PinkFoam, Monokote, Basswood, Balsa
from Aerothon.ACWingWeight import ACSolidWing, ACRibWing

#
# Create the wing
#
Wing = ACMainWing(1)
Wing.Lift_LO       = 60 * LBF    # Desired lift at liftoff
Wing.V_max_climb   = 65 * FT/SEC # liftoff speed
Wing.Alt_LO        = 200 * FT    # liftoff altitude
Wing.AR            = 8.0         # Aspect Ratio
Wing.b             = 12*12*IN       # Wingspan


###############################################################################
#
# Geometric properties
#
###############################################################################

Wing.FullWing = True

Wing.Fb      = [0.4,0.8,1] # Span fraction corresponding to inputs below
Wing.TR      = [1,0.8,.7]  # Taper ratio 
Wing.Gam     = [20*ARCDEG, 20*ARCDEG, 20*ARCDEG] # Dihedral
Wing.Lam     = [0*ARCDEG, 0*ARCDEG, 0*ARCDEG] # wing sweep angles
Wing.CEdge   = 'LE' # Defines constant edge
Wing.ConstUpper = True

###############################################################################
#
# Aerodynamic properties
#
###############################################################################

#
# Set the airfoils
#
Wing.Airfoil = 'CLi4005b' # Chooses airfoil from /AircraftDesign/Airfoils.. just load in normal .dat files like you would to xflr5
Wing.o_eff = 0.98 # Oswald Efficiency
Wing.FWCF = 0.98 # I'm... not srue what this does it isn't referenced in any calculations.

#
# Polar slope evaluations
#
Wing.ClSlopeAt = (6*ARCDEG, 7*ARCDEG) # Just the range that the polars will be drawn
Wing.CmSlopeAt = (-1*ARCDEG, 0*ARCDEG) # Range polars will be drawn

###############################################################################
#
# Control surfaces
#
###############################################################################

#
# Define the control surfaces
#
Wing.AddControl('Aileron')
Wing.Aileron.Fc = 0.25 # Percentage of chord aileron will take up
Wing.Aileron.Fb = 0.3 # Percentage of span
Wing.Aileron.Ft = 0.2 # Percentage distance from tip
Wing.Aileron.SgnDup = -1.

Wing.Aileron.Servo.Fc     = 0.3 # Pereentage of chord loction of servo
Wing.Aileron.Servo.Weight = 0.01*LBF

###############################################################################
#
# Structural properties
#
###############################################################################

#
# Spar material (basswood, 1/4in width at max airfoil thickness
#
sparw = 1.0*IN
spart = 1.0*IN
Basswood = Basswood.copy()
Balsa    = Balsa.copy()

BassLD = Basswood.ForceDensity * sparw * spart

#
# Rib material (1/8in balsa)
#
RibMat = Balsa.copy()
RibMat.Thickness = 0.125*IN

Wing.SetWeightCalc(ACRibWing)
Wing.WingWeight.AddSpar("MainSpar",1*IN,1*IN)
Wing.WingWeight.MainSpar.SparMat.LinearForceDensity = BassLD
Wing.WingWeight.SkinMat                             = Monokote.copy()
Wing.WingWeight.RibMat                              = RibMat
Wing.WingWeight.RibSpace                            = 5*IN

if __name__ == '__main__':
    import pylab as pyl
    
    print "V lift of   : ", AsUnit( Wing.GetV_LO(), 'ft/s' )
    print "V stall     : ", AsUnit( Wing.V_Stall, 'ft/s' )
    print "Wing Area   : ", AsUnit( Wing.S, 'in**2' )
    print "Wing Span   : ", AsUnit( Wing.b, 'ft' )
    print "Wing AR     : ", Wing.AR
    print "Wing MAC    : ", AsUnit( Wing.MAC(), 'in' )
    print "Wing Xac    : ", Wing.Xac()
    print "Wing dCM_da : ", Wing.dCM_da()
    print "Wing dCL_da : ", Wing.dCL_da()
    print "Lift of Load: ", AsUnit( Wing.Lift_LO, 'lbf' )

    print "Wing Thickness: ", Wing.Thickness(0*FT)
    print "Wing Chord    : ", Wing.Chord(0*FT)
    print "Wing Area     : ", Wing.S
    print "Wing Lift     : ", Wing.Lift_LO
    print
    print "Wing Weight : ", AsUnit( Wing.Weight, 'lbf' )
    print "Wing MOI    : ", AsUnit( Wing.MOI(), 'slug*ft**2' )
   
#    Wing.WriteAVLWing('MonoWing.avl')
    
   # Wing.Draw3DWingPolars(fig=3)
  #  Wing.Draw2DAirfoilPolars(fig=2)

    Wing.Draw(fig = 1)
    pyl.show()